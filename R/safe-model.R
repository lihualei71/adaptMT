#---------------------------------------------------------------
# Functions to fit models in adapt
#---------------------------------------------------------------

safe_glm <- function(formula, family, data, weights = NULL,
                      ...){
    options(warn = -1)

    formula <- as.formula(formula)
    if (family$link %in% c("inverse", "log")){
        fit <- try(glm(formula, family, data, weights, ...),
                   silent = TRUE)
        if (class(fit)[1] == "try-error"){
            mod_mat <- model.matrix(formula, data = data)
            p <- ncol(mod_mat) - 1
            start <- c(1, rep(0, p))
            fit <- glm(formula, family, data, weights,
                       start = start, ...)
        }
    } else {
        fit <- glm(formula, family, data, weights, ...)
    }

    fitv <- as.numeric(
        predict(fit, type = "response")
        )

    df <- fit$rank
    info <- list(df = df)

    options(warn = 0)
    
    return(list(fitv = fitv, info = info))
}

safe_gam <- function(formula, family, data, weights = NULL,
                      ...){
    options(warn = -1)
    
    formula <- as.formula(formula)
    if (family$link %in% c("inverse", "log")){
        fit <- try(mgcv::gam(formula, family, data, weights, ...),
                   silent = TRUE)
        if (class(fit)[1] == "try-error"){
            mod_mat <- model.matrix(formula, data = data)
            p <- ncol(mod_mat) - 1
            start <- c(1, rep(0, p))
            fit <- mgcv::gam(formula, family, data, weights,
                             start = start, ...)
        }
    } else {
        fit <- mgcv::gam(formula, family, data, weights, ...)
    }

    fitv <- as.numeric(
        predict(fit, type = "response")
        )

    df <- fit$rank
    info <- list(df = df)

    options(warn = 0)
    
    return(list(fitv = fitv, info = info))
}

safe_glmnet <- function(x, y, family, weights = NULL,
                        ...){
    options(warn = -1)

    if (class(family)[1] == "family"){
        family <- family$family
    }

    if (family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")){
        if (is.null(weights)){
            fit <- glmnet::cv.glmnet(x, y, 
                                     family = family, ...)
        } else {
            weights <- pminmax(weights, 1e-5, 1-1e-5)
            fit <- glmnet::cv.glmnet(x, y, weights,
                                     family = family, ...)
        }
    } else if (family == "Gamma"){
        if (is.null(weights)){
            fit <- HDtweedie::cv.HDtweedie(x, y, p = 2,
                                           standardize = TRUE,
                                           ...)
        } else {
            weights <- pminmax(weights, 1e-5, 1-1e-5)
            fit <- HDtweedie::cv.HDtweedie(x, y, p = 2,
                                           weights = weights,
                                           standardize = TRUE,
                                           ...)
        }
    }

    fitv <- as.numeric(
        predict(fit, newx = x, s = "lambda.min",
                type = "response")
        )

    beta <- coef(fit, s = "lambda.min")
    vi <- as.numeric(beta != 0)[-1]    
    df <- sum(vi) + 1
    info <- list(df = df, vi = vi)

    options(warn = 0)
    
    return(list(fitv = fitv, info = info))
}
