#---------------------------------------------------------------
# Compute oracle local FDR estimate by using revealing all
# p-values and compute the correlation between the oracle
# estimate and estimate within each step
#---------------------------------------------------------------

corr_lfdr <- function(obj, data = NULL){
    if (class(obj)[1] != "adapt"){
        stop("\'obj\; is not of class \'adapt\"")
    }
    
}
