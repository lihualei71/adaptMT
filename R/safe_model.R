#---------------------------------------------------------------
# Functions to fit models in adapt
#---------------------------------------------------------------

adapt_glm <- function(formula, family, data, weights = NULL,
                     ...){
    formula <- as.formula(formula)
    if (family$link %in% c("inverse", "log")){
        mod <- try(glm(formula, family, data, weights,
                       ...),
                   silent = TRUE)
        if (class(mod)[1] == "try-error"){
            mod_mat <- model.matrix(formula, data = data)
            p <- ncol(mod_mat) - 1
            start <- c(1, rep(0, p))
            mod <- glm(formula, family, data, weights,
                       start = start, ...)
        }
    } else if (family$link == "logit"){
        mod <- suppressWarnings(
            glm(formula, family, data, weights,
                ...))
    } else {
        mod <- glm(formula, family, data, weights,
                   ...)
    }

    fitv <- predict(mod, type = "response")

    return(list(mod = mod, fitv = fitv))
}

adapt_gam <- function(formula, family, data, weights = NULL,
                     ...){
    formula <- as.formula(formula)
    if (family$link %in% c("inverse", "log")){
        mod <- try(mgcv::gam(formula, family, data, weights,
                             ...),
                   silent = TRUE)
        if (class(mod)[1] == "try-error"){
            mod_mat <- model.matrix(formula, data = data)
            p <- ncol(mod_mat) - 1
            start <- c(1, rep(0, p))
            mod <- mgcv::gam(formula, family, data, weights,
                             start = start, ...)
        }
    } else if (family$link == "logit"){
        mod <- suppressWarnings(
            mgcv::gam(formula, family, data, weights,
                      ...))
    } else {
        mod <- mgcv::gam(formula, family, data, weights,
                         ...)
    }

    fitv <- predict(mod, type = "response")

    return(list(mod = mod, fitv = fitv))
}

adapt_glmnet <- function(x, y, family, weights = NULL,
                        ...){
    ## Safe GLMnet
    if (class(family)[1] == "family"){
        family <- family$family
    }

    if (family %in% c("gaussian", "poisson", "multinomial", "cox", "mgaussian")){
        mod <- glmnet::cv.glmnet(x, y, weights,
                                 family = family, ...)
    } else if (family == "binomial") {
        n <- length(y)
        newy <- c(rep(1, n), rep(0, n))
        newx <- rbind(x, x)
        weights <- c(y, 1 - y)
        mod <- glmnet::cv.glmnet(newx, newy, weights,
                                 family = "binomial", ...)
    } else if (family == "Gamma"){
        mod <- HDtweedie::cv.HDtweedie(x, y, p = 2,
                                       weights = weights,
                                       ...)
    }

    fitv <- as.numeric(
        predict(mod, newx = x, s = "lambda.min",
                type = "response")
        )
    return(list(mod = mod, fitv = fitv))
}
