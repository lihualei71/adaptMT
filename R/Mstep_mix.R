#===============================================================
# Compute M-step for the mixture model.
#===============================================================

Mstep_mix <- function(x, Hhat, phat, dist,
                      pifun = NULL, piargs = NULL,
                      mufun = NULL, muargs = NULL){
    input_type_pi <- func_input_type(pifun)
    input_type_mu <- func_input_type(mufun)
    response_name <- find_newname(colnames(x))

    if (input_type_pi== "formula"){
        data <- cbind(data.frame(Hhat), x)
        colnames(data)[1] <- response_name    
        piargs <- c(list(data = data), complete_formulas(piargs,
                             response_name))
    } else if (input_type_pi == "xy"){
        piargs <- list(x = x, y = Hhat)
    } else if (input_type_pi == "Xy"){
        piargs <- list(X = x, y = Hhat)
    } else {
        stop("pifun has irregular input types. Replace another function or writing a wrapper of pifun with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)")
    }

    if (input_type_mu== "formula"){
        data <- cbind(data.frame(dist$g(phat)), x)
        colnames(data)[1] <- response_name    
        muargs <- c(list(data = data), complete_formulas(muargs,
                             response_name))
    } else if (input_type_mu == "xy"){
        muargs <- list(x = x, y = Hhat)
    } else if (input_type_mu == "Xy"){
        muargs <- list(X = x, y = Hhat)
    } else {
        stop("mufun has irregular input types. Replace another function or writing a wrapper of mufun with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)")
    }

    fit_pi <- do.call(pifun, piargs)
    pix <- as.numeric(fit_pi$fitv)
    pix <- pmin(pmax(pix, 0), 1)
    
    fit_mu <- do.call(mufun, muargs)
    mux <- as.numeric(fit_mu$fitv)
    if (dist$family$family == "Gamma"){
        mux <- pmax(mux, 1)
    } else if (dist$family$family == "gaussian"){
        mux <- pmax(mux, 0)
    }

    if (any(is.na(pix))){
        stop("pix in M-step has NA values!")
    }
    if (any(is.na(mux))){
        stop("mux in M-step has NA values!")
    }
    return(list(pix = pix, mux = mux,
                fit_pi = fit_pi$mod, fit_mu = fit_mu$mod))
}

Mstep_mix_glm <- function(x, Hhat, phat, dist,
                          pi.formula, mu.formula,
                          input_type = "formula",
                          ...){
    if (!is.null(input_type) && (input_type != "formula")){
        input_type <- "formula"
        warning("\"input_type\" for GLM can only be formula!")
    }
    extra.args <- list(...)    
    pifun <- safe.logistic.glm
    piargs <- c(list(formula = pi.formula), extra.args)
    mufun <- function(formula, data, ...){
        safe.glm(formula, data,
                 family = dist$family, ...)
    }
    muargs <- c(list(formula = mu.formula), extra.args)
    res <- Mstep.mix(x = x, Hhat = Hhat, phat = phat, dist = dist,
                     input_type = input_type,
                     pifun = pifun,
                     piargs = piargs,
                     mufun = mufun,
                     muargs = muargs)
    pi.vi <- NA
    mu.vi <- NA

    pi.df <- res$fit_pi$rank
    mu.df <- res$fit_mu$rank
    df <- pi.df + mu.df

    res$fit_pi <- NULL
    res$fit_mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)
    return(res)    
}

Mstep_mix_gam <- function(x, Hhat, phat, dist,
                          pi.formula, mu.formula,
                          input_type = "formula",
                          ...){
    if (!is.null(input_type) && (input_type != "formula")){
        input_type <- "formula"
        warning("\"input_type\" for GAM can only be formula!")
    }
    extra.args <- list(...)    
    pifun <- safe.logistic.gam
    piargs <- c(list(formula = pi.formula), extra.args)
    mufun <- function(formula, data, ...){
        safe.gam(formula, data,
                 family = dist$family, ...)
    }
    muargs <- c(list(formula = mu.formula), extra.args)
    res <- Mstep.mix(x = x, Hhat = Hhat, phat = phat, dist = dist,
                     input_type = input_type,
                     pifun = pifun,
                     piargs = piargs,
                     mufun = mufun,
                     muargs = muargs)
    pi.vi <- NA
    mu.vi <- NA
        
    pi.df <- res$fit_pi$rank
    mu.df <- res$fit_mu$rank
    df <- pi.df + mu.df

    res$fit_pi <- NULL
    res$fit_mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)
    return(res)    
}

Mstep_mix_glmnet <- function(x, Hhat, phat, dist,
                             input_type = "xy",
                             ...){
    if (!is.null(input_type) && (input_type != "xy")){
        input_type <- "xy"
        warning("\"input_type\" for glmnet can only be xy!")
    }
    extra.args <- list(...)    
    pifun <- function(x, y, ...){
        y <- pmin(pmax(logit(y), logit(10^-15)), logit(1-10^-15))
        res <- safe.gaussian.glmnet(x = x, y = y, ...)
        res$fitv <- pmin(pmax(inv.logit(res$fitv), 10^-15), 1-10^-15)
        return(res)
    }
    piargs <- extra.args
    mufun <- function(x, y, ...){
        safe.glmnet(x, y,
                    family = dist$family, ...)
    }
    muargs <- extra.args
    res <- Mstep.mix(x = x, Hhat = Hhat, phat = phat, dist = dist,
                     input_type = input_type,
                     pifun = pifun,
                     piargs = piargs,
                     mufun = mufun,
                     muargs = muargs)
    pi.coef <- coef(res$fit_pi, s = "lambda.min")
    mu.coef <- coef(res$fit_mu, s = "lambda.min")

    pi.vi <- as.numeric(pi.coef != 0)[-1]
    mu.vi <- as.numeric(mu.coef != 0)[-1]
    
    pi.df <- sum(pi.vi)
    mu.df <- sum(mu.vi)
    df <- pi.df + mu.df

    res$fit_pi <- NULL
    res$fit_mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)    
    return(res)
}
