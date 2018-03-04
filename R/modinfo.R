modinfo <- function(fit){
    UseMethod("modinfo", fit)
}

modinfo.default <- function(fit){
    warning(paste0("No method is found for class(fit) = ", class(fit), " and output NA by default. Please specify modinfo.", class(fit), " to (at least) report a non-NA df (degree-of-freedom) if model selection is performed based on any information criterion, e.g. AIC, BIC. "))
    return(list(df = NA, vi = NA))
}

modinfo.glm <- modinfo.gam <- function(fit){
    df <- fit$rank
    vi <- NULL
    return(list(df = df, vi = vi))
}

modinfo.glmnet <- modinfo.HDtweedie <- function(fit){
    beta <- coef(fit$fit_pi, s = "lambda.min")
    vi <- as.numeric(beta != 0)[-1]    
    df <- sum(vi) + 1
    return(list(df = df, vi = vi))
}

