################################################################
## Compute the threshold curve s(x) given a level curve of local fdr. See Section 4.2.
##
## Required Input:
##    dist: distribution family for p-values in "exp_family" class.
##    params: parameters including pix and mux.
##    lfdr.lev: target local-fdr level.
## Output:
##    s: the threshold curve for a given local-fdr level.
################################################################

compute.threshold.mix <- function(dist, params, lfdr.lev){
    pix <- params$pix
    mux <- params$mux
    if (lfdr.lev == 0 || lfdr.lev == 1){
        return(rep(lfdr.lev, length(pix)))
    }
    
    val1 <- dist$f(1, mux) / lfdr.lev +
        (1 - pix) / pix * (1 - lfdr.lev) / lfdr.lev
    val2 <- (log(val1) + dist$A(mux) - dist$A(dist$mu.star)) /
        (dist$eta(mux) - dist$eta(dist$mu.star))
    dist$g.inv(val2)
}
