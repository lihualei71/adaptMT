#---------------------------------------------------------------
# EM algorithm to fit a mixture model.
#---------------------------------------------------------------

EM_loglik <- function(dist, params, Hhat, phat){
    pix <- params$pix
    mux <- params$mux
    loglik1 <- sum(Hhat * log(pix) + (1 - Hhat) * log(1 - pix))
    loglik2 <- sum(Hhat * log(dist$h(phat, mux)))
    return(loglik1 + loglik2)
}

EM_mix <- function(x, pvals, s, dist,
                   method = c("glm", "gam", "glmnet", "custom"),
                   algo = NULL, piargs = NULL, muargs = NULL,
                   params0 = list(pix = NULL, mux = NULL),
                   niter = 10, tol = 1e-4,
                   verbose = FALSE){
    method <- method[1]
    if (method == "custom" && is.null(algo)){
        stop("\"algo\" must be specified. The algo argument should be a list with names pifun, mufun, init_pifun, init_mufun")
    }

    if (verbose){
        cat("Model fitting starts!\n")
        cat("Reduce niter if it is too time-consuming.\n")
        pb <- utils::txtProgressBar(min = 0, max = niter, style = 3, width = 50)
    }

    if (is.null(params0$pix) || is.null(params0$mux)){
        init_res <- init_mix(x, pvals, s, dist, method,
                             algo$init_pifun, algo$init_mufun,
                             piargs, muargs)
        old_pix <- pix <- init_res$pix
        old_mux <- mux <- init_res$mux
    } else {
        old_pix <- pix <- params0$pix
        old_mux <- mux <- params0$mux
    }

    for (step in 1:niter){
        Estep_res <-
            Estep_mix(pvals, s, dist, pix, mux)
        Mstep_res <-
            Mstep_mix(x, Estep_res$Hhat, Estep_res$phat,
                      dist, method,
                      algo$pifun, algo$mufun,
                      piargs, muargs)
        pix <- Mstep_res$pix
        mux <- Mstep_res$mux
        if (max(abs(mux - old_mux)) < tol &&
            max(abs(pix - old_pix)) < tol){
            break
        }
        old_pix <- pix
        old_mux <- mux
        if (verbose){
            utils::setTxtProgressBar(pb, step)
        }        
    }
    if (verbose){    
        cat("\n")
    }
    params <- list(pix = pix, mux = mux)
    loglik <- EM_loglik(dist, params,
                        Estep_res$Hhat, Estep_res$phat)
    info <- Mstep_res$info
    
    return(list(params = params, loglik = loglik, info = info))
}
