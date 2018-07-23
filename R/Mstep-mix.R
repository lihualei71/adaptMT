#===============================================================
# Compute M-step for the mixture model.
#===============================================================

Mstep_mix <- function(x, pvals, dist,
                      Hhat, bhat,
                      pifun, mufun,
                      piargs = NULL, muargs = NULL,
                      type = "unweighted"){
    if (!"weights" %in% formalArgs(mufun)){
        stop("'mufun' does not have input 'weights'")
    }

    n <- length(Hhat)
    x_aug <- rbind(x, x)
    H_aug <- c(rep(1, n), rep(0, n))
    weights <- c(Hhat, 1 - Hhat)
    piargs <- complete_args(x_aug, H_aug, pifun, piargs, weights)
    pi_res <- fit_pi(pifun, piargs, type = "Mstep")
    pi_res$fitv <- pi_res$fitv[1:n]

    y_aug <- c(dist$g(pvals), dist$g(1 - pvals))

    if (type == "weighted"){
        weights <- c(Hhat * bhat, Hhat * (1 - bhat))
    } else if (type == "unweighted"){
        weights <- c(bhat, 1 - bhat)
    }
    muargs <- complete_args(x_aug, y_aug, mufun, muargs, weights)
    mu_res <- fit_mu(mufun, muargs, dist, type = "Mstep")
    mu_res$fitv <- mu_res$fitv[1:n]

    res <- list(pix = pi_res$fitv,
                mux = mu_res$fitv,
                pi_info = pi_res$info,
                mu_info = mu_res$info)

    return(res)
}
