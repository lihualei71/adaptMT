################################################################
## Compute E-step for the mixture model.
##
## Required Input:
##    x: covariate.
##    pvals: p-values.
##    plow: lower threshold, i.e. s(x) in AdaPT.
##    phigh: upper threshold, i.e. 1 - s(x) in AdaPT.
##    pix: pi(x) in last step.
##    mux: mu(x) in last step.
##    dist: distribution family for p-values in "exp_family" class.
## Output:
##    Hhat: conditional expectation of labels.
##    yhat: imputed y-values.
##    phat: g^{-1}(yhat).
################################################################

Estep.mix <- function(x, pvals, plow, phigh, dist,
                      pix, mux){
    fp <- dist$f(pvals, mux)
    fp.mirror <- dist$f(1 - pvals, mux)
    y <- dist$g(pvals)
    y.mirror <- dist$g(1 - pvals)
    Hhat <- ifelse(
        pvals < plow | pvals > phigh,
        1 / (1 + 2 * (1 - pix) / pix / (fp + fp.mirror)),
        1 / (1 + (1 - pix) / pix / fp))
    yhat <- ifelse(
        pvals < plow | pvals > phigh,
        (y * fp + y.mirror * fp.mirror) / (fp + fp.mirror),
        y)
    phat <- dist$g.inv(yhat)

    if (any(is.na(Hhat))){
        stop("Hhat in the E-step has NA values!")
    }
    if (any(is.na(phat))){
        stop("phat in the E-step has NA values!")
    }
    return(list(Hhat = Hhat, phat = phat))
}
