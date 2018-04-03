#===============================================================
# Compute E-step for the mixture model.
#===============================================================

## ' Computing E-step for Mixture Models
## '
## ' \code{Estep_mix} computes the E-step, namely imputation of missing values, including the p-values and indicators of null/non-null hypotheses.
## '
## ' The p-values are assumed to be generated from an covariate-varying mixture model
## ' \deqn{H_{i}\sim Ber(\pi(x_{i})), p_{i} \mid (H_{i} = 0) \sim U([0, 1]), p_{i} \mid (H_{i} = 1) \sim h(\cdot; \mu(x_{i}))}{Hi ~ Ber(\pi(xi)), pi | (Hi = 0) ~ U([0, 1]), pi | (Hi = 1) ~ h(p; \mu(xi))}
## ' where \eqn{h(p; \mu)} is the density of an exponential family (see \code{\link{exp_family}}). Given a threshold curve s(x), the partially i-th masked p-value is defined as \eqn{p_{i}}{pi} if \eqn{p_{i}\in [s(x_{i}), 1-s(x_{i})]}{pi in [s(xi), 1-s(xi)]}. The E-step computes the expectation of \eqn{H_{i}}{Hi} and \eqn{g(p_{i})H_{i}}{g(pi)Hi} given the parameters \eqn{\pi(x_{i})}{\pi(xi)} and \eqn{\mu(x_{i})}{\mu(xi)}, based on which computes the imputed values for \eqn{H_{i}}{Hi} and \eqn{p_{i}}{pi}. See ...
## '
## ' @param pvals a vector of values in [0, 1]. P-values
## ' @param s a vector of values in [0, 1]. Threshold curve
## ' @param dist an object of class "\code{\link{exp_family}}"
## ' @param pix a vector of values in [0, 1]. \eqn{\pi(x_{i})}{\pi(xi)}
## ' @param mux a vector of values. \eqn{\mu(x_{i})}{\mu(xi)}
## ' @return Imputed values, including
## ' \item{Hhat}{a vector of values in [0, 1]. Imputed values for \eqn{H_{i}}{Hi}'s}
## ' \item{phat}{a vector of values in [0, 1]. Imputed values for \eqn{p_{i}}{pi}'s}
Estep_mix <- function(pvals, s, dist, pix, mux){
    hp <- dist$h(pvals, mux)
    hp_mir <- dist$h(1 - pvals, mux)
    Hhat <- ifelse(
        pvals < s | pvals > 1 - s,
        1 / (1 + 2 * (1 - pix) / pix / (hp + hp_mir)),
        1 / (1 + (1 - pix) / pix / hp)
        )
    Hhat <- pminmax(Hhat, 1e-5, 1-1e-5)
    bhat <- ifelse(
        pvals < s | pvals > 1 - s,
        hp / (hp + hp_mir),
        1
        )

    if (any(is.na(Hhat))){
        stop("Hhat in the E-step has NAs.")
    }
    if (any(is.na(bhat))){
        stop("bhat in the E-step has NAs.")
    }

    return(list(Hhat = Hhat, bhat = bhat))
}
