#===============================================================
# Compute M-step for the mixture model.
#===============================================================

Mstep_mix_pi <- function(x, Hhat, fun, args){
    args <- complete_args(x, Hhat, fun, args)

    fit_pi(fun, args, type = "Mstep")
}

Mstep_mix_mu <- function(x, Hhat, phat, dist, fun, args){
    yhat <- dist$g(phat)
    Hhat <- pmax(Hhat, 1e-5)
    if (!"weights" %in% formalArgs(fun)){
        stop("mufun does not have input weights")
    }

    args <- complete_args(x, yhat, fun, args, Hhat)

    fit_mu(fun, args, dist, type = "Mstep")    
}

Mstep_mix <- function(x, Hhat, phat, dist,
                      pifun, mufun, 
                      piargs = NULL, muargs = NULL){
    pi_res <- Mstep_mix_pi(x, Hhat, pifun, piargs)
    mu_res <- Mstep_mix_mu(x, Hhat, phat, dist, mufun, muargs)
    c(pi_res, mu_res)
}
