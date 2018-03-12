#---------------------------------------------------------------
# AdaPT with varying-coefficient two-group model
# Lihua Lei & William Fithian,
#   "AdaPT: An interactive procedure for multiple testing with side information"
# Available from http://arxiv.org/abs/1609.06035
#---------------------------------------------------------------

fdp_hat <- function(A, R){
    (1 + A) / pmax(1, R)
}


compute_lfdr_mix <- function(pvals, dist, params){
    pix <- params$pix
    mux <- params$mux
    lfdr <- (pix * dist$h(1, mux) + 1 - pix) /
        (pix * dist$h(pvals, mux) + 1 - pix)
    return(lfdr)
}

compute_threshold_mix <- function(dist, params, lfdr_lev){
    pix <- params$pix
    mux <- params$mux
    if (lfdr_lev == 0 || lfdr_lev == 1){
        return(rep(lfdr_lev, length(pix)))
    }
    
    val1 <- dist$h(1, mux) / lfdr_lev +
        (1 - pix) / pix * (1 - lfdr_lev) / lfdr_lev
    val2 <- (log(val1) + dist$A(mux) - dist$A(dist$mustar)) /
        (dist$eta(mux) - dist$eta(dist$mustar))
    dist$ginv(val2)
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

#' @export
AdaPT <- function(x, pvals, models,
                  dist = beta_family(),
                  s0 = rep(0.45, length(pvals)),
                  alphas = seq(0.01, 1, 0.01),                  
                  params0 = list(pix = NULL, mux = NULL),
                  nfits = 20, nms = 1,
                  verbose = list(print = TRUE, fit = FALSE, ms = TRUE),
                  niter_fit = 10, tol = 1e-4,
                  niter_ms = 20, cr = "BIC",
                  return_data = TRUE){
    ## Check if "dist" is of class "exp_family"
    if (class(dist)[1] != "exp_family"){
        stop("\"dist\" must be of class \"exp_family\".")
    }

    ## When a single model is provided, set "nms" to be NULL
    if (class(models) == "adapt_model"){
        model <- models
        nms <- NULL
    } else if (is.list(models)){
        types <- sapply(models, class)
        if (any(types != "adapt_model")){
            stop("All models should be of class \"adapt_model\".")
        }
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
        niter = niter_fit, tol = tol,
        verbose = verbose$fit)
    if (any(stamps[, 2] == "ms")){
        ms_args_root <- list(
            x = x, pvals = pvals, dist = dist, models = models,
            cr = cr, niter = niter_ms, tol = tol,
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
    s_return <- matrix(0, n, m) # threshold
    params_return <- list() # parameters (including pix and mux)
    model_list <- list() # all selected models
    info_list <- list() # other information (df, vi, etc.)
    reveal_order <- which((pvals > s) & (pvals < 1 - s)) # the order to be revealed
    fdp_return <- rep(Inf, length(reveal_order)) # fdphat along the whole path

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
            model <- ms_res$model
            model_list <- append(model_list, list(model))
            info_list <- append(info_list, list(ms_res$info))
        } else if (type == "fit"){
            fit_args <- c(
                list(s = s, params0 = params, model = model),
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
        reveal_order <- c(reveal_order, inds)
        ## Shortcut to calculate FDPhat after revealing the hypotheses one by one
        Adecre <- cumsum(pvals[inds] >= 1 - s[inds])
        Rdecre <- cumsum(pvals[inds] <= s[inds])
        fdp <- fdp_hat(A - Adecre, R - Rdecre)
        fdp_return <- c(fdp_return, fdp)        
        fdp <- pmin(fdp, minfdp)
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
             s = s_return,
             reveal = data.frame(
                 order = reveal_order,
                 fdp = fdp_return
                 ),
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
