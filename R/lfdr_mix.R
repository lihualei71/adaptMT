################################################################
## Local fdr computation for the mixture model.
##
## Compute (over-estimated) Local fdr:
## lfdr_i = f(1|x)/f(p_i|x).
## 
## Required Input:
##    p: p-values.
##    dist: distribution family for p-values in "exp_family" class.
##    params: parameters including pix and mux.
##
## Output:
##    lfdr: local fdr for each p-value.
################################################################

lfdr.mix <- function(pvals, dist, params){
    pix <- params$pix
    mux <- params$mux
    lfdr <- (pix * dist$f(1, mux) + 1 - pix) /
        (pix * dist$f(pvals, mux) + 1 - pix)
    return(lfdr)
}
