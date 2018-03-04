#---------------------------------------------------------------
# AdaPT with varying-coefficient two-group model
# Lihua Lei & William Fithian,
#   "AdaPT: An interactive procedure for multiple testing with side information"
# Available from http://arxiv.org/abs/1609.06035
#---------------------------------------------------------------

fdp_hat <- function(A, R){
    (1 + A) / pmax(1, R)
}

AdaPT_preprocess <- function(methods, algos, piargs, muargs){
    ## Check the type of methods
    if (!class(methods) %in% c("NULL", "character", "list")){
        stop("\"methods\" must be NULL or a string or a list of strings.")
    } else if (class(methods) == "list" && !all(sapply(methods, is.character))){
        stop("\"methods\" must be a list of strings.") 
    }
    
    ## Check the type of algos
    if (!class(algos) %in% c("NULL", "list", "function")){
        stop("\"algos\" must be NULL or a function or a list of functions.")
    } else if (class(methods) == "list" && !all(sapply(methods, is.function))){
        stop("\"algos\" must be a list of functions.") 
    }

    ## Check the type of piargs/muargs. We allow piargs/muargs to be NULL, or a single argument list, or a list of argument lists.
    ## An argument list is a list with names/values being the argument names/values of the target function.
    ## When piargs/muargs is a list of argument lists, names(piargs)/names(muargs) should be NULL and
    ## each element should have the same set of argument names
    if (!is.null(piargs)){
        if (!is.list(piargs)){
            stop("\"piargs\" must be a list")
        } else if (!is.null(names(piargs))){
            if ("" %in% names(piargs)){
                stop("\"piargs\" has empty argument names.")
            } else {
                piargs <- list(piargs)
            }
        } else if (!all(sapply(piargs, is.list))){
            stop("\"piargs\" must be a list of (argument) lists")
        }
    }
    if (!is.null(muargs)){
        if (!is.list(muargs)){
            stop("\"muargs\" must be a list")
        } else if (!is.null(names(muargs))){
            if ("" %in% names(muargs)){
                stop("\"muargs\" has empty argument names.")
            } else {
                muargs <- list(muargs)
            }
        } else if (!all(sapply(muargs, is.list))){
            stop("\"muargs\" must be a list of (argument) lists")
        }
    }
    
    ## Check if the number of piargs equal to the number of muargs
    if (length(piargs) != length(muargs)){
        stop("\"piargs\" and \"muargs\" must have the same number of elements")
    }
    
    ## Check if either "methods" or "algos" is specified
    if (is.null(methods) && is.null(algos)){
        stop("Either \"methods\" or \"algos\" must be specified.")
    }
    
    ## If "methods" is unspecified but "algos" is specified, set it as a vector of "custom"s
    if (is.null(methods)){
        methods <- rep("custom", length(algos))
    }

    ## If "methods" is a list, transform it into a vector
    if (is.list(methods)){
        methods <- unlist(methods)
    }
    
    ## If "algos/piargs/muargs" is unspecified, set it as a list of lists
    if (is.null(algos)){
        algos <- lapply(1:length(methods), function(i) NULL)
    }
    if (is.null(piargs)){
        piargs <- lapply(1:length(methods), function(i) list())
    }
    if (is.null(muargs)){
        muargs <- lapply(1:length(methods), function(i) list())
    }

    ## When multiple methods are provided, check if the number matches that of piargs/muargs
    if (length(methods) > 1 && length(methods) != length(piargs)){
        stop("Number of methods, piargs and muargs should be equal.")
    }
    
    ## When a single method and a single set of piargs/muargs are provided, unlist algos/piargs/muargs
    if (length(methods) == 1 && length(piargs) == 1){
        algos <- algos[[1]]
        piargs <- piargs[[1]]
        muargs <- muargs[[1]]        
    }

    ## When a single method and multiple piargs/muargs are provided, set "methods" as a vector
    if (length(methods) == 1 && length(piargs) > 1){
        methods <- rep(methods, length(piargs))
    }

    return(list(methods = methods, algos = algos,
                piargs = piargs, muargs = muargs))
}

create_stamps <- function(nmasks, nfits, nms){
    fit_stamps <- c(seq(0, nmasks, floor(nmasks / nfits)))[1:nfits]
    stamps <- data.frame(stamp = fit_stamps, type = "fit",
                         stringsAsFactors = FALSE)
    if (!is.null(nms)){
        ms_inds <- seq(1, nfits, floor(nfits / nms))[1:nms]
        stamps[ms_inds, 2] <- "ms"
    }
    stamps <- rbind(stamps, data.frame(stamp = nmasks, type = "end"))
    return(stamps)
    
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

    ## Clean the inputs
    clean_inputs <- AdaPT_preprocess(methods, algos, piargs, muargs)
    methods <- clean_inputs$methods
    algos <- clean_inputs$algos
    piargs <- clean_inputs$piargs
    muargs <- clean_inputs$muargs
    
    ## When a single method is provided, wrap the arguments as "model" and set "nms" to be NULL
    if (length(methods) == 1){
        model <- list(method = methods, algo = algos,
                      piargs = piargs, muargs = muargs)
        nms <- NULL
    } else {
        if (!is.integer(nms) || nms <= 0){
            nms <- 1
        } else if (nms > nfits){
            warning("Model selection cannot be more frequent than model fitting. Set \"nms\" to \"nfits\"")
            nms <- nfits
        }
    }
    
    ## Create time stamps when model is fitted or model selection is performed
    nmasks <- sum(pvals <= s) + sum(pvals >= 1 - s)
    stamps <- create_stamps(nmasks, nfits, nms)    

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

        while (alphaind > 0){# check if lower FDR level is achieved
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
                        round(fdpnew, 4), ", Number of Rej. ",
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

        ## Update s(x)
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
