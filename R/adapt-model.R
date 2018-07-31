#===============================================================
# adapt_model class
#===============================================================

#' adapt_model Objects for M-steps
#'
#'
#' \code{adapt_model} objects provide the functions and their arguments in computing the M-steps.
#' Each object can be passed to \code{\link{adapt}} as a candidate model.
#'
#' \code{pifun} should be in the form of \code{pifun(formula, data, family, weights, ...)} or \code{pifun(x, y, family, ...)}.
#' The former includes \code{\link[stats]{glm}} and \code{\link[mgcv]{gam}} and the latter includes \code{\link[glmnet]{glmnet}}.
#' The outputs should be in the form of \code{list(fitv = , info = , ...)} where \code{fitv} gives the estimate of pi(x),
#' as a vector with the same order of \code{x}, and \code{info} should at least contain a key \code{df} if model selection is used, i.e. \code{info = list(df = , ...)}
#'
#' \code{mufun} should be in the form of \code{pifun(formula, data, family, weights, ...)} or \code{pifun(x, y, family, weights, ...)}.
#' Note that \code{mufun} must take \code{weights} as an input. The outputs should be in the same form as \code{pifun} except that \code{fitv} should give the estimate of mu(x).
#'
#' When \code{pifun} / \code{mufun} takes the form of \code{(formula, family, ...)}, \code{piargs} / \code{muargs} should at least contain a key \code{formula}; when \code{pifun} / \code{mufun} takes the form of \code{(x, y, family, ...)}, \code{piargs} / \code{muargs} can be empty.
#'
#' For glm/gam/glmnet, one can use the shortcut by running \code{\link{gen_adapt_model}} with name = "glm" or "gam" or "glmnet" but without specifying \code{pifun}, \code{mufun}, \code{pifun_init} and \code{mufun_init}. See examples below.
#'
#' @param pifun a function to fit pi(x). See Details
#' @param mufun a function to fit mu(x). See Details
#' @param pifun_init a function to fit pi(x) at the initial step
#' @param mufun_init a function to fit mu(x) at the initial step
#' @param piargs a list. Arguments for "pifun". An empty list as default
#' @param muargs a list. Arguments for "mufun". An empty list as default
#' @param piargs_init a list. Arguments for piargs_init. An empty list as default
#' @param muargs_init a list. Arguments for muargs_init. An empty list as default
#' @param name a string. An optional argument for the user-specified name of the model. An empty string as default.
#'
#' @return
#' \item{name}{same as the input \code{name}}
#' \item{algo}{a list recording \code{pifun}, \code{mufun}, \code{pifun_init} and \code{mufun_init}}
#' \item{args}{a list recording \code{piargs}, \code{muargs}, \code{piargs_init} and \code{muargs_init}}
#'
#' @examples
#' \donttest{
#' # Exemplary code to generate 'adapt_model' for logistic-Gamma glm  with naive initialization.
#' # The real implementation in the package is much more complicated.
#'
#' # pifun as a logistic regression
#' pifun <- function(formula, data, weights, ...){
#'   glm(formula, data, weights = weights, family = binomial(),  ...)
#' }
#' # pifun_init as a constant
#' pifun_init <- function(x, pvals, s, ...){
#'   rep(0.1, length(pvals))
#' }
#' # mufun as a Gamma GLM
#' mufun <- function(formula, data, weights, ...){
#'   glm(formula, data, weights = weights, family = Gamma(), ...)
#' }
#' # mufun_init as a constant
#' mufun_init <- function(x, pvals, s, ...){
#'   rep(1.5, length(pvals))
#' }
#'
#' library("splines") # for using ns() in the formula
#' piargs <- list(formula = "ns(x, df = 8)")
#' muargs <- list(formula = "ns(x, df = 8)")
#' name <- "glm"
#'
#' mod <- gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
#'                        piargs, muargs, name = name)
#' mod
#'
#' # Using shortcut for GLM. See the last paragraph of Details.
#' mod2 <- gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' mod2
#' }
#'
#' @export
gen_adapt_model <- function(pifun = NULL,
                            mufun = NULL,
                            pifun_init = NULL,
                            mufun_init = NULL,
                            piargs = list(),
                            muargs = list(),
                            piargs_init = list(),
                            muargs_init = list(),
                            name = ""){
    args <- list(piargs = piargs, muargs = muargs,
                 piargs_init = piargs_init, muargs_init = muargs_init)

    if (is.null(pifun) && is.null(mufun) &&
        is.null(pifun_init) && is.null(mufun_init)){
        model <- structure(
            list(name = name, args = args),
            class = "adapt_model"
            )
        return(model)
    }

    if (!is.function(pifun)){
        stop("\"pifun\" must be a function.")
    }
    if (!is.function(mufun)){
        stop("\"mufun\" must be a function.")
    }
    if (!is.function(pifun_init)){
        stop("\"pifun_init\" must be a function.")
    }
    if (!is.function(mufun_init)){
        stop("\"mufun_init\" must be a function.")
    }

    algo <- list(pifun = pifun, mufun = mufun,
                 pifun_init = pifun_init, mufun_init = mufun_init)

    model <- structure(
        list(name = name, algo = algo, args = args),
        class = "adapt_model"
        )
    return(model)
}

gen_adapt_model_glm <- function(dist,
                                piargs = list(),
                                muargs = list()){
    pifun <- function(formula, data, weights, ...){
        safe_glm(formula, data, weights = weights,
                 family = quasibinomial(), ...)
    }

    mufun <- function(formula, data, weights, ...){
        safe_glm(formula, data, weights = weights,
                 family = dist$family, ...)
    }

    pifun_init <- function(x, pvals, s, ...){
        J <- ifelse(
            pvals < s | pvals > 1 - s, 1,
            2 * s / (2 * s - 1)
            )
        if (any(s >= 0.49)){
            J[s >= 0.49] <- 0
        }

        fun <- function(formula, data, ...){
            safe_glm(formula, data, family = gaussian(), ...)
        }

        args <- list(...)
        args <- complete_args(x, J, fun, args)

        fit_pi(fun, args, type = "init")
    }

    mufun_init <- function(x, pvals, s, ...){
        phat <- ifelse(
            pvals < s | pvals > 1 - s,
            pmin(pvals, 1 - pvals),
            pvals
            )
        phat <- pminmax(phat, 1e-15, 1-1e-15)
        yhat <- dist$g(phat)

        fun <- function(formula, data, ...){
            safe_glm(formula, data, family = dist$family, ...)
        }

        args <- list(...)
        args <- complete_args(x, yhat, fun, args)

        fit_mu(fun, args, dist, type = "init")
    }

    if (is.null(piargs$formula) || is.null(muargs$formula)){
        stop("Argument \"formula\" is missing from \"piargs\" or \"muargs\".")
    }

    piargs_init <- piargs
    muargs_init <- muargs

    gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
                    piargs, muargs, piargs_init, muargs_init,
                    name = "glm")
}

gen_adapt_model_gam <- function(dist,
                                piargs = list(),
                                muargs = list()){
    pifun <- function(formula, data, weights, ...){
        safe_gam(formula, data, weights = weights,
                 family = quasibinomial(), ...)
    }

    mufun <- function(formula, data, weights, ...){
        safe_gam(formula, data, weights = weights,
                 family = dist$family, ...)
    }

    pifun_init <- function(x, pvals, s, ...){
        J <- ifelse(
            pvals < s | pvals > 1 - s, 1,
            2 * s / (2 * s - 1)
            )
        if (any(s >= 0.49)){
            J[s >= 0.49] <- 0
        }

        fun <- function(formula, data, ...){
            safe_gam(formula, data, family = gaussian(), ...)
        }

        args <- list(...)
        args <- complete_args(x, J, fun, args)

        fit_pi(fun, args, type = "init")
    }

    mufun_init <- function(x, pvals, s, ...){
        phat <- ifelse(
            pvals < s | pvals > 1 - s,
            pmin(pvals, 1 - pvals),
            pvals
            )
        phat <- pminmax(phat, 1e-15, 1-1e-15)
        yhat <- dist$g(phat)

        fun <- function(formula, data, ...){
            safe_gam(formula, data, family = dist$family, ...)
        }

        args <- list(...)
        args <- complete_args(x, yhat, fun, args)

        fit_mu(fun, args, dist, type = "init")
    }

    if (is.null(piargs$formula) || is.null(muargs$formula)){
        stop("Argument \"formula\" is missing from \"piargs\" or \"muargs\".")
    }

    piargs_init <- piargs
    muargs_init <- muargs

    gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
                    piargs, muargs, piargs_init, muargs_init,
                    name = "gam")
}

gen_adapt_model_glmnet <- function(dist,
                                   piargs = list(),
                                   muargs = list()){
    pifun <- function(x, y, weights, ...){
        safe_glmnet(x, y, weights = weights,
                    family = "binomial", ...)
    }


    mufun <- function(x, y, weights, ...){
        safe_glmnet(x, y, weights = weights,
                    family = dist$family, ...)
    }

    pifun_init <- function(x, pvals, s, ...){
        J <- ifelse(
            pvals < s | pvals > 1 - s, 1,
            2 * s / (2 * s - 1)
            )
        if (any(s >= 0.49)){
            J[s >= 0.49] <- 0
        }

        fun <- function(x, y, ...){
            safe_glmnet(x, y, family = "gaussian", ...)
        }

        args <- list(...)
        args <- complete_args(x, J, fun, args)

        fit_pi(fun, args, type = "init")
    }

    mufun_init <- function(x, pvals, s, ...){
        phat <- ifelse(
            pvals < s | pvals > 1 - s,
            pmin(pvals, 1 - pvals),
            pvals
            )
        phat <- pminmax(phat, 1e-15, 1-1e-15)
        yhat <- dist$g(phat)

        fun <- function(x, y, ...){
            safe_glmnet(x, y, family = dist$family$family, ...)
        }

        args <- list(...)
        args <- complete_args(x, yhat, fun, args)

        fit_mu(fun, args, dist, type = "init")
    }

    piargs_init <- piargs
    muargs_init <- muargs

    gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
                    piargs, muargs, piargs_init, muargs_init,
                    name = "glmnet")
}

