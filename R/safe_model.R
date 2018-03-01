################################################################
## Safe versions of model fitting.
## Maximally circumvent unexpected errors during the process.
##
## Required Input:
##     data: a data-frame
## Output:
##     mod: model
##     fitv: fitted value
##
## Main modification:
## Add an argument "alter.formulas" to incorporate a list of
## alternative models when the fitting fails. For instance,
## the default model to use is a GLM with natural spline basis
## ns(x, df = 4) and the alternative models are GLMs with natural
## spline basis ns(x, df = 3), ns(x, df = 2), ns(x, df = 1). This
## effectively avoids the potential colinearity problem during the
## process.
################################################################

#' Securing model-fitting algorithm (internal function)
#'
#' A generic handle for making a fitting algorithm secure to run by adding a list of
#' candidate arguments. When the fitting function encounters an error, it prints the
#' current error message and proceeds to run the function on the candidate argument,
#' one-by-one from the list. It stops when a particular set of arguments does not give
#' an error and the results of which will be outputed. If all sets of candidate arguments
#' report errors, the function will report an error message, in addition to the previous
#' one.
#'
#' @param algo A function (fitting algorithm).
#' @param data_args A list, with the keys being the name of an argument of \code{algo} for specifying
#'     data/input. Examples include list(data = ...) or list(x = ..., y = ...).
#' @param algo_args A list, with the keys being the name of an argument of \code{algo}, other than
#'     those specified for data/input, and the values being the value of that argument. Examples
#'     include list(formula = ...).
#' @param alter_args A list of list, with each sub-list of the same lengths and keys as
#'     \code{algo_args}, but with different values. NULL if no candidate argument.
#' @return One of the following: \itemize{
#'    \item The output of \code{algo}, when running \code{algo} on \code{algo_args} does not report an error,
#'    \item an error message indicating the error by running \code{algo} on \code{algo_args}, and the output
#'    of \code{algo} on the first entry of \code{alter_args} which does not report an error,
#'    \item an error message indicating the error by running \code{algo} on \code{algo_args}, together with
#'    an error message indicating that none of candidate arguments work.
#'    }
safe_fit <- function(algo, data_args, algo_args,
                     alter_args = NULL,
                     ...){
    old <- options(warn = -1)
    on.exit(options(old), add = TRUE)

    extra_args <- list(...)
    args <- c(data_args, algo_args, extra_args)
    mod <- try(do.call(algo, args))

    if (class(mod)[1] != "try-error" || is.null(alter_args)){
        return(list(mod = mod, args = args))
    }

    rm(args)
    gc()
    print(mod$condition)

    for (algo_args in alter_args){
        args <- c(data_args, algo_args, extra_args)
        mod <- try(do.call(algo, args))

        if (class(mod)[1] != "try-error"){
            return(list(mod = mod, args = args))
        }

        rm(args)
        gc()
    }

    stop("model fitting fails!")
}

safe_glm <- function(formula, data,
                     alter.formulas = NULL,
                     family = gaussian(),
                     ...){
    ## Safe GLM
    glm <- stats::glm
    algo <- function(formula, data, family, ...){
        if (family$link %in% c("inverse", "log")){
            mod <- try(glm(formula = formula, data = data,
                           family = family, ...),
                       silent = TRUE)
            if (class(mod)[1] != "try-error"){
                return(mod)
            }
            tmp.mat <- model.matrix(formula, data = data)
            p <- ncol(tmp.mat) - 1
            start <- c(1, rep(0, p))
            mod <- glm(formula = formula, data = data,
                       family = family, start = start, ...)
        } else {
            mod <- glm(formula = formula, data = data,
                       family = family, ...)
        }
        return(mod)
    }

    data_args <- c(list(data = data))
    algo_args <- c(list(formula = formula, family = family))
    alter_args <- lapply(alter.formulas, function(formula){
        c(list(formula = formula, family = family))
    })

    result <- safe_fit(algo, data_args, algo_args, alter_args,
                       ...)
    mod <- result$mod
    fitv <- predict(mod, type = "response")

    return(list(mod = mod, fitv = fitv))
}

safe_gam <- function(formula, data,
                     alter.formulas = NULL,
                     family = gaussian(),
                     ...){
    algo <- function(formula, data, family, ...){
        if (family$link %in% c("inverse", "log")){
            mod <- try(mgcv::gam(formula = formula, data = data,
                           family = family, ...),
                       silent = TRUE)
            if (class(mod)[1] != "try-error"){
                return(mod)
            }
            tmp.mat <- model.matrix(formula, data = data)
            p <- ncol(tmp.mat) - 1
            start <- c(1, rep(0, p))
            mod <- mgcv::gam(formula = formula, data = data,
                             family = family, start = start, ...)
        } else {
            mod <- mgcv::gam(formula = formula, data = data,
                             family = family, ...)
        }
        return(mod)
    }

    data_args <- c(list(data = data))
    algo_args <- c(list(formula = formula, family = family))
    alter_args <- lapply(alter.formulas, function(formula){
        c(list(formula = formula, family = family))
    })

    result <- safe_fit(algo, data_args, algo_args, alter_args,
                       ...)
    mod <- result$mod
    fitv <- predict(mod, type = "response")

    return(list(mod = mod, fitv = fitv))
}

safe_glmnet <- function(x, y,
                        family = gaussian(),
                        ...){
    ## Safe GLMnet
    if (class(family)[1] == "family"){
        family <- family$family
    }

    if (family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")){
        algo <- function(...)
            glmnetUtils::cv.glmnet(..., family = family)
    } else if (family == "Gamma"){
        algo <- function(...)
            HDtweedie::cv.HDtweedie(..., p = 2)
    }

    data_args <- list(x = x, y = y)
    algo_args <- list()
    alter_args <- NULL
    result <- safe_fit(algo, data_args, algo_args, alter_args, ...)
    mod <- result$mod
    fitv <- as.numeric(predict(mod, newx = x, s = "lambda.min", type = "response"))
    return(list(mod = mod, fitv = fitv))
}

safe_logistic_glm <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe.glm(formula, data, alter.formulas,
             family = binomial(),
             ...)
}

safe_logistic_gam <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe_gam(formula, data, alter.formulas,
             family = binomial(),
             ...)
}

safe_logistic_glmnet <- function(x, y,
                                 ...){
    safe_glmnet(x, y,
                family = "binomial",
                ...)
}

safe_gaussian_glm <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe_glm(formula, data, alter.formulas,
             family = gaussian(),
             ...)
}

safe_gaussian_gam <- function(formula, data,
                              alter.formulas = NULL,
                              ...){
    safe_gam(formula, data, alter.formulas,
             family = gaussian(),
             ...)
}

safe_gaussian_glmnet <- function(x, y,
                                 ...){
    safe_glmnet(x, y,
                family = "gaussian",
                ...)
}

safe_gamma_glm <- function(formula, data,
                           alter.formulas = NULL,
                           ...){
    safe_glm(formula, data, alter.formulas,
             family = Gamma(),
             ...)
}

safe_gamma_gam <- function(formula, data,
                           alter.formulas = NULL,
                           ...){
    safe_gam(formula, data, alter.formulas,
             family = Gamma(),
             ...)
}

safe_gamma_glmnet <- function(x, y,
                              ...){
    safe_glmnet(x, y,
                family = "Gamma",
                ...)
}
