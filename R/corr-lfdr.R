#---------------------------------------------------------------
# Compute oracle local FDR estimate by revealing all
# p-values and compute the correlation between the oracle
# estimate and estimate within each step
#---------------------------------------------------------------

#' Quantifying Information Loss of Adaptive P-Value Thresholding
#'
#' \code{corr_lfdr} computes the oracle local FDR estimate, by using revealing all p-values, and computes the Pearson correlation between it and the estimate within each step of \code{adapt}.
#'
#' @param obj an 'adapt' object. Output of \code{\link{adapt}} function
#' @param x covariates (i.e. side-information). Should be compatible to \code{models}.
#' @param pvals a vector of values in [0, 1]. P-values
#' @param model an optional argument. If \code{model = NULL} then the last model in \code{obj$models} is used for fitting the oracle model (i.e. with all p-values revealed). Otherwise it should be an 'adapt_model' object
#' @param niter_oracle an positive integer. Number of iterations in EM algorithm
#'
#' @return
#' \itemize{
#' \item{corr}{a vector of values in [0, 1]. Pearson correlation of oracle local FDR estimate and the estimates within each step. Each value corresponds to an entry of \code{obj$params}}
#' \item{oracle_lfdr}{a vector of values in [0, 1]. Oracle local FDR estimate}
#' \item{lfdr}{a matrix of values in [0, 1]. Local FDR estimates within each step.}
#' \item{alphas}{a vector of values in [0, 1]. The target FDR levels corresponding to each local FDR estimate}
#' \item{nmasks}{a vector of integers. The number of masked p-values corresponding to each local FDR estimate}
#' }
#'
#' @examples
#' \donttest{
#' # Load estrogen data
#' data(estrogen)
#' pvals <- as.numeric(estrogen$pvals)
#' x <- data.frame(x = as.numeric(estrogen$ord_high))
#' dist <- beta_family()
#'
#' # Subsample the data for convenience
#' inds <- (x$x <= 5000)
#' pvals <- pvals[inds]
#' x <- x[inds,,drop = FALSE]
#'
#' # Run adapt_glm
#' library("splines")
#' formulas <- paste0("ns(x, df = ", 6:10, ")")
#' res <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
#'                  mu_formulas = formulas, dist = dist, nfits = 10)
#'
#' # Run corr_lfdr
#' obj <- corr_lfdr(res, x, pvals)
#' obj$corr
#' }
#'
#' @export
corr_lfdr <- function(obj, x, pvals, model = NULL,
                      niter_oracle = 100){
    if (class(obj)[1] != "adapt"){
        stop("\'obj\' is not of class \'adapt\"")
    }
    n <- length(pvals)

    params_list <- obj$params
    dist <- obj$dist
    m <- length(params_list)
    alphas <- rep(0, m)
    nmasks <- rep(0, m)
    lfdr <- matrix(0, n, m)
    for (i in 1:m){
        params <- list(pix = params_list[[i]]$pix,
                    mux = params_list[[i]]$mux)
        alphas[i] <- params_list[[i]]$alpha
        nmasks[i] <- params_list[[i]]$nmasks
        lfdr[, i] <- compute_lfdr_mix(pvals, dist, params)
    }

    ## Oracle lfdr
    if (is.null(model)){
        model <- obj$models[[length(obj$models)]]
    }
    oracle_params <- EM_mix(x, pvals, rep(0, n), dist, model,
                            params0 = params,
                            niter = niter_oracle)$params
    oracle_lfdr <- compute_lfdr_mix(pvals, dist, oracle_params)

    ## Correlation
    corr <- apply(lfdr, 2, function(x){cor(x, oracle_lfdr)})

    return(list(corr = corr,
                oracle_lfdr = oracle_lfdr,
                lfdr = lfdr,
                alphas = alphas, nmasks = nmasks))
}
