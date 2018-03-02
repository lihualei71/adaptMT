#---------------------------------------------------------------
# Local fdr computation for the mixture model.
#---------------------------------------------------------------

lfdr_mix <- function(pvals, dist, params){
    pix <- params$pix
    mux <- params$mux
    lfdr <- (pix * dist$h(1, mux) + 1 - pix) /
        (pix * dist$h(pvals, mux) + 1 - pix)
    return(lfdr)
}
