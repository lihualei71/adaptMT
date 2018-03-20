#---------------------------------------------------------------
# Functions to fit models in adapt
#---------------------------------------------------------------

safe_glm <- function(formula, family, data, weights = NULL,
                      ...){
    oldw <- options(warn = -1)

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
    } else if (family$link == "logit"){
        fit <- glm(formula, family, data, weights, ...)
    } else {
        fit <- glm(formula, family, data, weights, ...)
    }

    fitv <- as.numeric(
        predict(fit, type = "response")
        )

    df <- fit$rank
    info <- list(df = df)

    options(oldw)
    
    return(list(fitv = fitv, info = info))
}

safe_gam <- function(formula, family, data, weights = NULL,
                      ...){
    oldw <- options(warn = -1)
    
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
    } else if (family$link == "logit"){
        fit <- mgcv::gam(formula, family, data, weights, ...)
    } else {
        fit <- mgcv::gam(formula, family, data, weights, ...)
    }

    fitv <- as.numeric(
        predict(fit, type = "response")
        )

    df <- fit$rank
    info <- list(df = df)

    options(oldw)
    
    return(list(fitv = fitv, info = info))
}

safe_glmnet <- function(x, y, family, weights = NULL,
                         ...){
    oldw <- options(warn = -1)
    
    if (class(family)[1] == "family"){
        family <- family$family
    }

    if (family %in% c("gaussian", "poisson", "multinomial", "cox", "mgaussian")){
        fit <- glmnet::cv.glmnet(x, y, weights,
                                 family = family, ...)
    } else if (family == "binomial") {
        n <- length(y)
        newy <- c(rep(1, n), rep(0, n))
        newx <- rbind(x, x)
        weights <- c(y, 1 - y)
        fit <- glmnet::cv.glmnet(newx, newy, weights,
                                 family = "binomial", ...)
    } else if (family == "Gamma"){
        fit <- HDtweedie::cv.HDtweedie(x, y, p = 2,
                                       weights = weights,
                                       ...)
    }

    fitv <- as.numeric(
        predict(fit, newx = x, s = "lambda.min",
                type = "response")
        )

    beta <- coef(fit, s = "lambda.min")
    vi <- as.numeric(beta != 0)[-1]    
    df <- sum(vi) + 1
    info <- list(df = df, vi = vi)

    options(oldw)
    
    return(list(fitv = fitv, info = info))
}
