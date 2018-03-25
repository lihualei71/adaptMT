#===============================================================
# Compute M-step for the mixture model.
#===============================================================

Mstep_mix <- function(x, pvals, dist,
                      Hhat, bhat, 
                      pifun, mufun, 
                      piargs = NULL, muargs = NULL){
    if (!"weights" %in% formalArgs(mufun)){
        stop("'mufun' does not have input 'weights'")
    }

    n <- length(Hhat)
    x_aug <- rbind(x, x)
    H_aug <- c(rep(1, n), rep(0, n))
    weights <- c(Hhat, 1 - Hhat) 
    piargs <- complete_args(x_aug, H_aug, pifun, piargs, weights)
    pi_res <- fit_pi(pifun, piargs, type = "Mstep")
    pi_res$pix <- pi_res$pix[1:n]

    y_aug <- c(dist$g(pvals), dist$g(1 - pvals))
    weights <- c(Hhat * bhat, Hhat * (1 - bhat))
    muargs <- complete_args(x_aug, y_aug, mufun, muargs, weights)
    mu_res <- fit_mu(mufun, muargs, dist, type = "Mstep")
    mu_res$mux <- mu_res$mux[1:n]

    return(c(pi_res, mu_res))
}
