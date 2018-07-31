#---------------------------------------------------------------
# Wrappers of AdaPT
#---------------------------------------------------------------

check_formulas <- function(formulas){
    err <- FALSE
    if (is.character(formulas)){
        formulas <- as.list(formulas)
    } else if (is.list(formulas)){
        err <- any(!sapply(formulas, class) %in% c("character", "formula"))
    } else {
        err <- TRUE
    }

    if (err){
        stop("Invalid formulas")
    }

    return(formulas)
}

#' Adaptive P-value Thresholding with Generalized Linear Models
#'
#' \code{adapt_glm} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[stats]{glm}}.
#'
#' \code{pi_formulas} and \code{mu_formulas} can either be a list or a vector with each element being a string or a formula. For instance, suppose \code{x} has a single column with name \code{x1}, the following five options are valid for the same inputs (\code{\link[splines]{ns}} forms a spline basis with \code{df} knots):
#' \enumerate{
#' \item{c("x1", "ns(x1, df = 8)");}
#' \item{c("~ x1", "~ ns(x1, df = 8)");}
#' \item{list("x1", "ns(x1, df = 8)");}
#' \item{list("~ x1", "~ ns(x1, df = 8)");}
#' \item{list(~ x1, ~ ns(x1, df = 8))}
#' }
#' There is no need to specify the name of the response variable, as this is handled in the function.
#'
#' When \code{x} has a few variables, it is common to use non-parametric GLM by replacing \code{x} by a spline basis of \code{x}. In this case, \code{\link[splines]{ns}} from \code{library(splines)} package is suggested.
#' 
#' @param pi_formulas a vector/list of strings/formulas. Formulas for fitting pi(x) by glm. See Details
#' @param mu_formulas a vector/list of strings/formulas. Formulas for fitting mu(x) by glm. See Details
#' @param piargs a list. Other arguments passed to glm for fitting pi(x)
#' @param muargs a list. Other arguments passed to glm for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}} (except \code{models})
#' @inheritParams adapt
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
#' # Run adapt by manually setting models for glm
#' models <- lapply(formulas, function(formula){
#'     piargs <- muargs <- list(formula = formula)
#'     gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' })
#' res2 <- adapt(x = x, pvals = pvals, models = models,
#'               dist = dist, nfits = 10)
#' 
#' # Check equivalence
#' identical(res, res2)
#' }
#'
#' @seealso
#' \code{\link{adapt}}, \code{\link{adapt_gam}}, \code{\link{adapt_glmnet}}, \code{\link[stats]{glm}}, \code{\link[splines]{ns}}
#' 
#' @export
adapt_glm <- function(x, pvals, pi_formulas, mu_formulas,
                      dist = beta_family(),
                      s0 = rep(0.45, length(pvals)),
                      alphas = seq(0.01, 1, 0.01),
                      piargs = list(), muargs = list(),
                      ...){
    if (!is.data.frame(x)){
        stop("\'x\' must be a data.frame")
    }

    pi_formulas <- check_formulas(pi_formulas)
    mu_formulas <- check_formulas(mu_formulas)
    stopifnot(length(pi_formulas) == length(mu_formulas))
        
    models <- lapply(1:length(pi_formulas), function(i){
        piargs <- c(list(formula = pi_formulas[[i]]), piargs)
        muargs <- c(list(formula = mu_formulas[[i]]), muargs)
        gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
    })

    adapt(x, pvals, models, dist, s0, alphas, ...)
}

#' Adaptive P-value Thresholding with Generalized Additive Models
#'
#' \code{adapt_gam} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[mgcv]{gam}} from \code{mgcv} package.
#'
#' \code{pi_formulas} and \code{mu_formulas} can either be a list or a vector with each element being a string or a formula. For instance, suppose \code{x} has a single column with name \code{x1}, the following five options are valid for the same inputs (\code{\link[splines]{ns}} forms a spline basis with \code{df} knots and \code{\link[mgcv]{s}} forms a spline basis with knots automatically selected by generalized cross-validation):
#' \enumerate{
#' \item{c("x1", "ns(x1, df = 8)", "s(x1)");}
#' \item{c("~ x1", "~ ns(x1, df = 8)", "s(x1)");}
#' \item{list("x1", "ns(x1, df = 8)", "s(x1)");}
#' \item{list("~ x1", "~ ns(x1, df = 8)", "s(x1)");}
#' \item{list(~ x1, ~ ns(x1, df = 8), s(x1))}
#' }
#' There is no need to specify the name of the response variable, as this is handled in the function.
#'
#' When \code{x} has a few variables, it is common to use non-parametric GLM by replacing \code{x} by a spline basis of \code{x}. In this case, \code{\link[splines]{ns}} from \code{library(splines)} package or \code{\link[mgcv]{s}} from \code{mgcv} package are suggested. When \code{\link[mgcv]{s}} (from \code{mgcv} package) is used, it is treated as a single model because the knots will be selected automatically.
#' 
#' @param pi_formulas a vector/list of strings/formulas. Formulas for fitting pi(x) by gam. See Details
#' @param mu_formulas a vector/list of strings/formulas. Formulas for fitting mu(x) by gam. See Details
#' @param piargs a list. Other arguments passed to gam for fitting pi(x)
#' @param muargs a list. Other arguments passed to gam for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}} (except \code{models})
#' @inheritParams adapt
#'
#' @seealso
#' \code{\link{adapt}}, \code{\link{adapt_glm}}, \code{\link{adapt_glmnet}}, \code{\link[mgcv]{gam}}, \code{\link[splines]{ns}}, \code{\link[mgcv]{s}}
#'
#' @examples
#' \donttest{
#' # Generate a 2-dim x
#' n <- 400
#' x1 <- x2 <- seq(-100, 100, length.out = 20)
#' x <- expand.grid(x1, x2)
#' colnames(x) <- c("x1", "x2")
#'
#' # Generate p-values (one-sided z test)
#' # Set all hypotheses in the central circle with radius 30 to be
#' # non-nulls. For non-nulls, z~N(2,1) and for nulls, z~N(0,1).
#' H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
#' mu <- ifelse(H0, 2, 0)
#' set.seed(0)
#' zvals <- rnorm(n) + mu
#' pvals <- 1 - pnorm(zvals)
#'
#' # Run adapt_gam with a 2d spline basis
#' library("mgcv")
#' formula <- "s(x1, x2)"
#' dist <- beta_family()
#' res <- adapt_gam(x = x, pvals = pvals, pi_formulas = formula,
#'                  mu_formulas = formula, dist = dist, nfits = 5)
#' }
#' 
#' 
#' @export
adapt_gam <- function(x, pvals, pi_formulas, mu_formulas,
                      piargs = list(), muargs = list(),
                      dist = beta_family(),
                      s0 = rep(0.45, length(pvals)),
                      alphas = seq(0.01, 1, 0.01),
                      ...){
    if (!is.data.frame(x)){
        stop("\'x\' must be a data.frame")
    }
    
    if (!requireNamespace("mgcv", quietly = TRUE)){
        stop("'mgcv' package is required for 'adapt_gam'. Please intall.")
    }

    pi_formulas <- check_formulas(pi_formulas)
    mu_formulas <- check_formulas(mu_formulas)
    stopifnot(length(pi_formulas) == length(mu_formulas))
    
    models <- lapply(1:length(pi_formulas), function(i){
        piargs <- c(list(formula = pi_formulas[[i]]), piargs)
        muargs <- c(list(formula = mu_formulas[[i]]), muargs)
        gen_adapt_model(name = "gam", piargs = piargs, muargs = muargs)
    })

    adapt(x, pvals, models, dist, s0, alphas, ...)
}

#' Adaptive P-value Thresholding with L1/L2 Penalized Generalized Linear Models
#'
#' \code{adapt_glmnet} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[glmnet]{glmnet}} from \code{glmnet} package.
#'
#' \code{adapt_glmnet} by default implements LASSO on \code{x} with lambda selected by cross-validation. Specify in \code{piargs} and \code{muargs} if ridge or elastic-net penalty is needed.
#' 
#' @param piargs a list. Other arguments passed to glmnet for fitting pi(x)
#' @param muargs a list. Other arguments passed to glmnet for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}} (except \code{models})
#' @inheritParams adapt
#'
#' 
#' @seealso
#' \code{\link{adapt}}, \code{\link{adapt_glm}}, \code{\link{adapt_gam}}, \code{\link[glmnet]{glmnet}}
#'
#' @examples
#' \donttest{
#' # Generate a 100-dim covariate x
#' set.seed(0)
#' m <- 100
#' n <- 1000
#' x <- matrix(runif(n * m), n, m)
#'
#' # Generate the parameters from a conditional two-group
#' # logistic-Gamma GLM  where pi(x) and mu(x) are both
#' # linear in x. pi(x) has an intercept so that the average
#' # of pi(x) is 0.3
#' inv_logit <- function(x) {exp(x) / (1 + exp(x))}
#' pi1 <- 0.3
#' beta.pi <- c(3, 3, rep(0, m-2))
#' beta0.pi <- uniroot(function(b){
#'     mean(inv_logit(x %*% beta.pi + b)) - pi1
#' }, c(-100, 100))$root
#' pi <- inv_logit(x %*% beta.pi + beta0.pi)
#' beta.mu <- c(2, 2, rep(0, m-2))
#' beta0.mu <- 0
#' mu <- pmax(1, x %*% beta.mu + beta0.mu)
#'
#' # Generate p-values
#' H0 <- as.logical(ifelse(runif(n) < pi, 1, 0))
#' y <- ifelse(H0, rexp(n, 1/mu), rexp(n, 1))
#' pvals <- exp(-y)
#'
#' # Run adapt_glmnet
#' res <- adapt_glmnet(x, pvals, s0 = rep(0.15, n), nfits = 5)
#' }
#' @export
adapt_glmnet <- function(x, pvals,
                         piargs = list(), muargs = list(),
                         dist = beta_family(),
                         s0 = rep(0.45, length(pvals)),
                         alphas = seq(0.01, 1, 0.01),
                         ...){
    if (!is.matrix(x) && !inherits(x, "sparseMatrix")){
        stop("Invalid \'x\'. See \'?glmnet\' for details.")
    }
    
    if (!requireNamespace("glmnet", quietly = TRUE)){
        stop("'glmnet' package is required for 'adapt_glmnet'. Please intall.")
    }

    models <- gen_adapt_model(name = "glmnet", piargs = piargs, muargs = muargs)

    adapt(x, pvals, models, dist, s0, alphas, ...)
}
