#---------------------------------------------------------------
# EM algorithm to fit a mixture model.
#---------------------------------------------------------------

EM_loglik <- function(pvals, dist, pix, mux, Hhat, bhat){
    loglik1 <- sum(Hhat * log(pix) + (1 - Hhat) * log(1 - pix))
    loglik2 <- sum(Hhat * bhat * log(dist$h(pvals, mux)) +
                   Hhat * (1 - bhat) * log(dist$h(1 - pvals, mux)))
    return(loglik1 + loglik2)
}

EM_mix <- function(x, pvals, s, dist, model,
                   params0 = list(pix = NULL, mux = NULL),
                   niter = 10, tol = 1e-4,
                   verbose = FALSE,
                   type = "unweighted"){
    model <- complete_model(model, dist)
    if (verbose){
        cat("Model fitting starts!\n")
        cat("Reduce niter if it is too time-consuming.\n")
        pb <- utils::txtProgressBar(min = 0, max = niter, style = 3, width = 50)
    }

    if (is.null(params0$pix) || is.null(params0$mux)){
        piargs_init <- c(list(x = x, pvals = pvals, s = s),
                         model$args$piargs_init)
        pix <- do.call(model$algo$pifun_init, piargs_init)$fitv
        muargs_init <- c(list(x = x, pvals = pvals, s = s),
                         model$args$muargs_init)
        mux <- do.call(model$algo$mufun_init, muargs_init)$fitv

            ## init_res <- init_mix(
            ## x, pvals, s, dist,
            ## model$algo$pifun_init, model$algo$mufun_init,
            ## model$args$piargs_init, model$args$muargs_init)

        old_pix <- pix ## <- init_res$pix
        old_mux <- mux ## <- init_res$mux
    } else {
        old_pix <- pix <- params0$pix
        old_mux <- mux <- params0$mux
    }

    for (step in 1:niter){
        Estep_res <-
            Estep_mix(pvals, s, dist, pix, mux)
        Mstep_res <-
            Mstep_mix(x, pvals, dist,
                      Estep_res$Hhat, Estep_res$bhat, 
                      model$algo$pifun, model$algo$mufun,
                      model$args$piargs, model$args$muargs,
                      type = type[1])
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
    loglik <- EM_loglik(pvals, dist, params$pix, params$mux,
                        Estep_res$Hhat, Estep_res$bhat)
    info <- list(pi = Mstep_res$pi_info, mu = Mstep_res$mu_info)
    
    return(list(params = params, loglik = loglik, info = info))
}
