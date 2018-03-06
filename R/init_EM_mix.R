#===============================================================
# Initialization of EM
#===============================================================

init_mix_pi <- function(x, pvals, s, fun, args){
    J <- ifelse(
        pvals < s | pvals > 1 - s, 1,
        2 * s / (2 * s - 1)
        )
    if (any(s >= 0.49)){
        J[s >= 0.49] <- 0
    }

    args <- complete_args(x, J, fun, args)

    fit_pi(fun, args, type = "init")
}

init_mix_mu <- function(x, pvals, s, dist, fun, args){
    phat <- ifelse(
        pvals < s | pvals > 1 - s,
        pmin(pvals, 1 - pvals),
        pvals
        )
    phat <- pminmax(phat, 1e-15, 1-1e-15)
    yhat <- dist$g(phat)

    args <- complete_args(x, yhat, fun, args)

    fit_mu(fun, args, dist, type = "init")
}

init_mix <- function(x, pvals, s, dist,
                     pifun, mufun,
                     piargs = NULL, muargs = NULL){
    pi_res <- init_mix_pi(x, pvals, s, pifun, piargs)
    mu_res <- init_mix_mu(x, pvals, s, dist, mufun, muargs)
    c(pi_res, mu_res)
}
