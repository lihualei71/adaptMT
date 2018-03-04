#---------------------------------------------------------------
# AdaPT with varying-coefficient two-group model
# Lihua Lei & William Fithian,
#   "AdaPT: An interactive procedure for multiple testing with side information"
# Available from http://arxiv.org/abs/1609.06035
#---------------------------------------------------------------

fdp_hat <- function(A, R){
    (1 + A) / pmax(1, R)
}

AdaPT <- function(x, pvals,
                  dist = beta_family(),
                  methods = NULL, algos = NULL,
                  piargs = NULL, muargs = NULL,
                  s0 = rep(0.45, length(pvals)),
                  alphas = seq(0.01, 1, 0.01),                  
                  params0 = list(pix = NULL, mux = NULL),
                  nfits = 20, nms = 1,
                  verbose = list(print = TRUE, fit = FALSE, ms = TRUE),
                  num_steps_fit = 10, tol = 1e-4,
                  num_steps_ms = 20, cr = "BIC",
                  return_data = TRUE){
    ## Check if "dist" is of class "exp_family"
    if (class(dist)[1] != "exp_family"){
        stop("\"dist\" must be of class \"exp_family\".")
    }

    ## Check if either "methods" or "algos" is specified
    if (is.null(methods) && is.null(algos)){
        stop("Either \"methods\" or \"algos\" must be specified.")
    }

    ## If "methods" is unspecified, set it as a vector of "custom"s
    ## If "methods" is list, transform it into a vector
    if (is.null(methods)){
        methods <- rep("custom", length(algos))
    } else if (is.list(methods)){
        methods <- unlist(methods)
    }

    ## If "algos/piargs/muargs" is unspecified, set it as a list of lists
    if (is.null(algos)){
        algos <- lapply(1:length(methods), function(i) list())
    }
    if (is.null(piargs)){
        piargs <- lapply(1:length(methods), function(i) list())
    }
    if (is.null(muargs)){
        muargs <- lapply(1:length(methods), function(i) list())
    }

    ## When a single method is provided, transform algo into a function (instead of a list of functions)
    if (length(methods) == 1 && is.list(algos[[1]])){
        algos <- algos[[1]]
    }

    ## When a single method and multiple piargs/muargs is provided, transform algo into a function (instead of a list of functions)
    if (length(methods) == 1 && is.list(piargs[[1]])){
        piargs <- piargs[[1]]
    }
    if (length(methods) == 1 && is.list(muargs[[1]])){
        muargs <- muargs[[1]]
    }

    ## When multiple methods are provided, check if the number matches that of piargs/muargs
    lens <- c(length(methods), length(piargs), length(muargs))
    if (length(methods) > 1 && max(lens) > min(lens)){
        stop("Number of methods, piargs and muargs should be equal.")
    }
    
    ## When a single method is provided, wrap the arguments "methods", "algos", "piargs" and "muargs" as "model"
    if (length(methods) == 1){
        model <- list(method = methods, algo = algos,
                      piargs = piargs, muargs = muargs)
    }
    
    ## Create time stamps when model is fitted or model selection is performed
    nmasks <- sum(pvals <= s) + sum(pvals >= 1 - s)
    fit_stamps <- c(seq(0, nmasks, floor(nmasks / nfits)))[1:nfits]
    stamps <- data.frame(stamp = fit_stamps, type = "fit",
                         stringsAsFactors = FALSE)
    if (!is.null(nms) && length(methods) > 1){
        ms_inds <- seq(1, nfits, floor(nfits / nms))[1:nms]
        stamps[ms_inds, 2] <- "ms"
    }
    stamps <- rbind(stamps, data.frame(stamp = nmasks, type = "end"))
    

    ## Create root arguments to simplify fitting and model selection
    fit_args_root <- list(
        x = x, pvals = pvals, dist = dist,
        num_steps = num_steps_fit, tol = tol,
        verbose = verbose$fit)
    if (any(stamps[, 2] == "ms")){
        ms_args_root <- list(
            x = x, pvals = pvals, dist = dist, cr = cr,
            cand_methods = methods,
            cand_algos = algos,
            cand_piargs = piargs,
            cand_muargs = muargs,
            num_steps = num_steps_ms, tol = tol,
            verbose = verbose$ms
            )
    }
    
    ## Initialization
    n <- length(pvals)
    params <- params0
    s <- s0
    A <- sum(pvals >= 1 - s)
    R <- sum(pvals <= s)
    minfdp <- fdp_hat(A, R) # initial FDPhat

    ## Remove the alphas greater than the initial FDPhat, except the smallest one among them
    alphas <- sort(alphas)
    if (min(alphas) >= minfdp){
        warning("Task completed! Initial \"s0\" has guaranteed FDR control for all alpha's in \"alphas\".")
        m <- 1
        alphaind <- 0
            
    } else if (max(alphas) < minfdp){
        m <- length(alphas)
        alphaind <- m
    } else {
        m <- max(which(alphas <= minfdp)) + 1
        alphaind <- m - 1
    }
    alphas <- alphas[1:m]

    nrejs_return <- rep(0, m) # number of rejections
    fdp_return <- rep(0, m) # post-hoc fdp estimate for each alpha
    s_return <- matrix(0, n, m) # threshold
    params_return <- list() # parameters (including pix and mux)
    model_list <- list() # all selected models
    info_list <- list() # other information (df, vi, etc.)
    qvals <- rep(NA, n) # q-values
    
    if (m > alphaind){
        nrejs_return[m] <- R
        fdp_return[m] <- minfdp
        s_return[, m] <- s0        
    }

    for (i in 1:(nrow(stamps) - 1)){
        if (alphaind == 0){
            if (verbose$print){
                cat("Task completed!")
            }
            break
        }

        start <- stamps[i, 1]
        end <- stamps[i + 1, 1]
        nreveals <- end - start # number of hypotheses to be revealed
        type <- stamps[i, 2] # type of model update
        mask <- (pvals <= s) | (pvals >= 1 - s)
        A <- sum(pvals >= 1 - s)
        R <- sum(pvals <= s)
        
        ## Model selection or model fitting
        if (type == "ms"){
            ms_args <- c(
                list(s = s, params0 = params),
                ms_args_root
                )
            ## Use "EM_mix_ms" from "EM-mix-ms.R"
            ms_res <- do.call(EM_mix_ms, ms_args)
            params <- ms_res$params
            model <- ms_res$best
            model_list <- append(model_list, list(model))
            info_list <- append(info_list, list(ms_res$info))
        } else if (type == "fit"){
            fit_args <- c(
                list(s = s, params0 = params,
                     method = model$method,
                     algo = model$algo,
                     piargs = model$piargs,
                     muargs = model$muargs),
                fit_args_root
                )
            ## Use "EM_mix" from "EM-mix.R"            
            fit_res <- do.call(EM_mix, fit_args)
            params <- fit_res$params
            model_list <- append(model_list, model)
            info_list <- append(info_list, list(fit_res$info))
        }
        params_return <- append(params_return, params)

        ## Estimate local FDR
        lfdr <- compute_lfdr_mix(
            pmin(pvals, 1 - pvals),
            dist, params)
        ## Find the top "nreveals" hypotheses with highest lfdr 
        lfdr[!mask] <- -Inf
        inds <- order(lfdr, decreasing = TRUE)[1:nreveals]
        ## Shortcut to calculate FDPhat after revealing the hypotheses one by one
        Adecre <- cumsum(pvals[inds] >= 1 - s[inds])
        Rdecre <- cumsum(pvals[inds] <= s[inds])
        fdp <- fdp_hat(A - Adecre, R - Rdecre)
        fdp <- pmin(fdp, minfdp)
        ## Calculate q-values
        qvals[inds] <- cummin(fdp)
        ## Calculate the current minimum FDPhat
        minfdp <- min(fdp)

        while (alphaind > 0){
            alpha <- alphas[alphaind]
            if (any(fdp <= alpha)){
                breakpoint <- which(fdp <= alpha)[1]
                lfdr_lev <- lfdr[inds[breakpoint]]
                snew <- compute_threshold_mix(dist, params, lfdr_lev)
                snew <- pmin(s, snew)
                s_return[, alphaind] <- snew
                
                fdpnew <- fdp[breakpoint]
                Rnew <- sum(pvals <= snew)
                fdp_return[alphaind] <- fdpnew
                nrejs_return[alphaind] <- Rnew

                if (verbose$print){
                    cat(paste0(
                        "alpha = ", alpha, ": FDPhat ",
                        round(fdpnew, 3), ", Number of Rej. ",
                        Rnew, "\n"))
                }

                alphaind <- alphaind - 1
            } else {
                break
            }
        }
        
        if (alphaind == 0){ # check again to save computation
            if (verbose$print){
                cat("Task completed!")
            }
            break
        }

        final_lfdr_lev <- lfdr[tail(inds, 1)]
        snew <- compute_threshold_mix(dist, params, final_lfdr_lev)
        s <- pmin(s, snew)
    }

    res <- structure(
        list(nrejs = nrejs_return,
             fdp = fdp_return,
             s = s_return,
             qvals = qvals,
             params = params_return,
             alphas = alphas,
             model = model_list,
             info = info_list),
        class = "AdaPT")
    if (return_data){
        res$data <- list(x = x, pvals = pvals)
    }
    return(res)
}

AdaPT.glm <- function(x, pvals, dist,
                      pi.formulas, mu.formulas,
                      glm.args = NULL,
                      num.steps.up = 10,
                      num.steps.ms = 30,
                      tol = 1e-4,
                      ...){
    if (class(pi.formulas)[1] == "list"){
        pi.formulas <- unlist(pi.formulas)
    }
    if (class(mu.formulas)[1] == "list"){
        mu.formulas <- unlist(mu.formulas)
    }
    nf <- length(pi.formulas)
    if (nf != length(mu.formulas)){
        stop("the numbers of formulas for pi(x) and mu(x) do not match!")
    }
    if (nf == 1){
        gen.algo.args <- c(list(pi.formula = pi.formulas,
                                mu.formula = mu.formulas),
                           glm.args)
        algo <- do.call(gen.AdaPT.glm, gen.algo.args)
        cand_algos <- NULL
    } else {
        algo <- NULL
        cand_algos <- lapply(1:nf, function(i){
            gen.algo.args <- c(list(pi.formula = pi.formulas[i],
                                    mu.formula = mu.formulas[i]),
                               glm.args)
            do.call(gen.AdaPT.glm, gen.algo.args)
        })
    }
    fit_args <- list(num.steps = num.steps.up, tol = tol)
    ms_args <- list(num.steps = num.steps.ms, tol = tol)
    ms.fun <- EM.mix.ms
    AdaPT(x = x, pvals = pvals, dist = dist,
          algo = algo, fit_args = fit_args,
          ms.fun = ms.fun, cand_algos = cand_algos,
          ms_args = ms_args, ...)
}

AdaPT.gam <- function(x, pvals, dist,
                      pi.formulas, mu.formulas,
                      gam.args = NULL,
                      num.steps.up = 10,
                      num.steps.ms = 30,
                      tol = 1e-4,
                      ...){
    if (class(pi.formulas)[1] == "list"){
        pi.formulas <- unlist(pi.formulas)
    }
    if (class(mu.formulas)[1] == "list"){
        mu.formulas <- unlist(mu.formulas)
    }
    nf <- length(pi.formulas)
    if (nf != length(mu.formulas)){
        stop("Error: the numbers of formulas for pi(x) and mu(x) do not match!")
    }
    if (nf == 1){
        gen.algo.args <- c(list(pi.formula = pi.formulas,
                                mu.formula = mu.formulas),
                           gam.args)
        algo <- do.call(gen.AdaPT.gam, gen.algo.args)
        cand_algos <- NULL
    } else {
        algo <- NULL
        cand_algos <- lapply(1:m, function(i){
            gen.algo.args <- c(list(pi.formula = pi.formulas[i],
                                    mu.formula = mu.formulas[i]),
                               gam.args)
            do.call(gen.AdaPT.gam, gen.algo.args)
        })
    }
    fit_args <- list(num.steps = num.steps.up, tol = tol)
    ms_args <- list(num.steps = num.steps.ms, tol = tol)
    ms.fun <- EM.mix.ms
    AdaPT(x = x, pvals = pvals, dist = dist,
          algo = algo, fit_args = fit_args,
          ms.fun = ms.fun, cand_algos = cand_algos,
          ms_args = ms_args, ...)
}

AdaPT.glmnet <- function(x, pvals, dist,
                         glmnet.args = NULL,
                         num.steps = 10,
                         tol = 1e-4,
                         ...){
    algo <- do.call(gen.AdaPT.glmnet, glmnet.args)
    cand_algos <- NULL
    fit_args <- list(num.steps = num.steps, tol = tol)
    ms.fun <- EM.mix.ms
    AdaPT(x = x, pvals = pvals, dist = dist,
          algo = algo, fit_args = fit_args,
          ms.fun = ms.fun, cand_algos = cand_algos, ...)
}

