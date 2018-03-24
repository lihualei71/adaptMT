#' Adaptive P-value Thresholding with Generalized Linear Models
#'
#' \code{adapt_glm} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[stats]{glm}}.
#'
#' \code{pi_formulas} and \code{mu_formulas}
#' 
#' @param pi_formulas a vector/list of strings/formulas. Formulas for fitting pi(x) by glm. See Details
#' @param mu_formulas a vector/list of strings/formulas. Formulas for fitting mu(x) by glm. See Details
#' @param piargs a list. Other arguments passed to glm for fitting pi(x)
#' @param muargs a list. Other arguments passed to glm for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}}
#' @inheritParams adapt
#'
#' @examples
#' \dontrun{
#' # Load estrogen data
#' data(estrogen)
#' pvals <- as.numeric(estrogen$pvals)
#' x <- data.frame(x = as.numeric(estrogen$ord))
#' dist <- beta_family()
#'
#' # Subsample the data for convenience
#' pvals <- pvals[1:5000]
#' x <- x[1:5000,,drop = FALSE]
#' 
#' # Run adapt_glm
#' library("splines")
#' formulas <- paste0("ns(x, df = ", 6:10, ")")
#' res <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
#'                  mu_formulas = formulas, dist = dist)
#'
#' # Run adapt by manually setting models for glm
#' models <- lapply(formulas, function(formula){
#'     piargs <- muargs <- list(formula = formula)
#'     gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' })
#' res2 <- adapt(x = x, pvals = pvals, models = models, dist = dist)
#' 
#' # Check equivalence
#' identical(res, res2)
#' }
#' 
#' @export
adapt_glm <- function(x, pvals, pi_formulas, mu_formulas,
                      dist = beta_family(),
                      s0 = rep(0.45, length(pvals)),
                      alphas = seq(0.01, 1, 0.01),
                      piargs = list(), muargs = list(),
                      ...){
    models <- lapply(1:length(pi_formulas), function(i){
        piargs <- c(list(formula = pi_formulas[[i]]), piargs)
        muargs <- c(list(formula = mu_formulas[[i]]), muargs)
        gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
    })

    adapt(x, pvals, models, dist, s0, alphas, ...)
}

#' Adaptive P-value Thresholding with Generalized Additive Models
#' 
#' @param pi_formulas a vector/list of strings/formulas. Formulas for fitting pi(x) by gam. See Details
#' @param mu_formulas a vector/list of strings/formulas. Formulas for fitting mu(x) by gam. See Details
#' @param piargs a list. Other arguments passed to gam for fitting pi(x)
#' @param muargs a list. Other arguments passed to gam for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}}
#' @inheritParams adapt
#' 
#' @export
adapt_gam <- function(x, pvals, pi_formulas, mu_formulas,
                      piargs = list(), muargs = list(),
                      dist = beta_family(),
                      s0 = rep(0.45, length(pvals)),
                      alphas = seq(0.01, 1, 0.01),
                      ...){
    if (!requireNamespace("mgcv", quietly = TRUE)){
        stop("'mgcv' package is required for 'adapt_gam'. Please intall.")
    }
    models <- lapply(1:length(pi_formulas), function(i){
        piargs <- c(list(formula = pi_formulas[[i]]), piargs)
        muargs <- c(list(formula = mu_formulas[[i]]), muargs)
        gen_adapt_model(name = "gam", piargs = piargs, muargs = muargs)
    })

    adapt(x, pvals, models, dist, s0, alphas, ...)
}

#' Adaptive P-value Thresholding with L1/L2 Penalized Generalized Linear Models
#' 
#' @param piargs a list. Other arguments passed to glmnet for fitting pi(x)
#' @param muargs a list. Other arguments passed to glmnet for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}}
#' @inheritParams adapt
#' 
#' @export
adapt_glmnet <- function(x, pvals,
                         piargs = list(), muargs = list(),
                         dist = beta_family(),
                         s0 = rep(0.45, length(pvals)),
                         alphas = seq(0.01, 1, 0.01),
                         ...){
    if (!requireNamespace("glmnet", quietly = TRUE)){
        stop("'glmnet' package is required for 'adapt_glmnet'. Please intall.")
    }

    models <- gen_adapt_model(name = "glmnet", piargs = piargs, muargs = muargs)

    adapt(x, pvals, models, dist, s0, alphas, ...)
}
