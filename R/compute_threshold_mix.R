#---------------------------------------------------------------
# Compute thresholds s(x) given a level curve of local fdr. 
#---------------------------------------------------------------

compute_threshold_mix <- function(dist, params, lfdr_lev){
    pix <- params$pix
    mux <- params$mux
    if (lfdr_lev == 0 || lfdr_lev == 1){
        return(rep(lfdr_lev, length(pix)))
    }
    
    val1 <- dist$h(1, mux) / lfdr_lev +
        (1 - pix) / pix * (1 - lfdr_lev) / lfdr_lev
    val2 <- (log(val1) + dist$A(mux) - dist$A(dist$mustar)) /
        (dist$eta(mux) - dist$eta(dist$mustar))
    dist$ginv(val2)
}
