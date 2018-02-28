################################################################
## Mstep for the mixture model.
## 
## Required Input:
##    x: covariate.
##    Hhat, phat: output of the E-step.
##    dist: distribution family for p-values in "exp_family" class.
##    pi.fit.fun: a function to fit initial pix using \tilde{J_i}.
##    pi.fit.args: other arguments passed into pi.fit.fun.
##    mu.fit.fun: a model to fit initial mux using (x,min(p,1-p)).
##    mu.fit.args: other arguments passed into mu.fit.fun.       
## Output:
##    pix, mux
################################################################

Mstep.mix <- function(x, Hhat, phat, dist,
                      input.type = c("formula", "xy"),
                      pi.fit.fun = NULL,
                      pi.fit.args = NULL,
                      mu.fit.fun = NULL,
                      mu.fit.args = NULL){
    input.type <- input.type[1]

    if (input.type == "formula"){
        response.name <- find.newname(colnames(x))
        data <- cbind(data.frame(Hhat), x)
        colnames(data)[1] <- response.name    
        pi.fit.args <- c(list(data = data),
                         complete.formulas(pi.fit.args,
                                           response.name))

        data <- cbind(data.frame(dist$g(phat)), x)
        colnames(data)[1] <- response.name

        mu.fit.args <- c(list(data = data),
                         complete.formulas(mu.fit.args,
                                           response.name))
    } else if (input.type == "xy"){
        pi.fit.args <- list(x = x, y = Hhat)
        mu.fit.args <- list(x = x, y = dist$g(phat))
    }

    fit.pi <- do.call(pi.fit.fun, pi.fit.args)
    pix <- as.numeric(fit.pi$fitv)
    pix <- pmin(pmax(pix, 0), 1)
    
    fit.mu <- do.call(mu.fit.fun, mu.fit.args)
    mux <- as.numeric(fit.mu$fitv)
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
                fit.pi = fit.pi$mod, fit.mu = fit.mu$mod))
}

Mstep.mix.glm <- function(x, Hhat, phat, dist,
                          pi.formula, mu.formula,
                          input.type = "formula",
                          ...){
    if (!is.null(input.type) && (input.type != "formula")){
        input.type <- "formula"
        warning("\"input.type\" for GLM can only be formula!")
    }
    extra.args <- list(...)    
    pi.fit.fun <- safe.logistic.glm
    pi.fit.args <- c(list(formula = pi.formula), extra.args)
    mu.fit.fun <- function(formula, data, ...){
        safe.glm(formula, data,
                 family = dist$family, ...)
    }
    mu.fit.args <- c(list(formula = mu.formula), extra.args)
    res <- Mstep.mix(x = x, Hhat = Hhat, phat = phat, dist = dist,
                     input.type = input.type,
                     pi.fit.fun = pi.fit.fun,
                     pi.fit.args = pi.fit.args,
                     mu.fit.fun = mu.fit.fun,
                     mu.fit.args = mu.fit.args)
    pi.vi <- NA
    mu.vi <- NA

    pi.df <- res$fit.pi$rank
    mu.df <- res$fit.mu$rank
    df <- pi.df + mu.df

    res$fit.pi <- NULL
    res$fit.mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)
    return(res)    
}

Mstep.mix.gam <- function(x, Hhat, phat, dist,
                          pi.formula, mu.formula,
                          input.type = "formula",
                          ...){
    if (!is.null(input.type) && (input.type != "formula")){
        input.type <- "formula"
        warning("\"input.type\" for GAM can only be formula!")
    }
    extra.args <- list(...)    
    pi.fit.fun <- safe.logistic.gam
    pi.fit.args <- c(list(formula = pi.formula), extra.args)
    mu.fit.fun <- function(formula, data, ...){
        safe.gam(formula, data,
                 family = dist$family, ...)
    }
    mu.fit.args <- c(list(formula = mu.formula), extra.args)
    res <- Mstep.mix(x = x, Hhat = Hhat, phat = phat, dist = dist,
                     input.type = input.type,
                     pi.fit.fun = pi.fit.fun,
                     pi.fit.args = pi.fit.args,
                     mu.fit.fun = mu.fit.fun,
                     mu.fit.args = mu.fit.args)
    pi.vi <- NA
    mu.vi <- NA
        
    pi.df <- res$fit.pi$rank
    mu.df <- res$fit.mu$rank
    df <- pi.df + mu.df

    res$fit.pi <- NULL
    res$fit.mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)
    return(res)    
}

Mstep.mix.glmnet <- function(x, Hhat, phat, dist,
                             input.type = "xy",
                             ...){
    if (!is.null(input.type) && (input.type != "xy")){
        input.type <- "xy"
        warning("\"input.type\" for glmnet can only be xy!")
    }
    extra.args <- list(...)    
    pi.fit.fun <- function(x, y, ...){
        y <- pmin(pmax(logit(y), logit(10^-15)), logit(1-10^-15))
        res <- safe.gaussian.glmnet(x = x, y = y, ...)
        res$fitv <- pmin(pmax(inv.logit(res$fitv), 10^-15), 1-10^-15)
        return(res)
    }
    pi.fit.args <- extra.args
    mu.fit.fun <- function(x, y, ...){
        safe.glmnet(x, y,
                    family = dist$family, ...)
    }
    mu.fit.args <- extra.args
    res <- Mstep.mix(x = x, Hhat = Hhat, phat = phat, dist = dist,
                     input.type = input.type,
                     pi.fit.fun = pi.fit.fun,
                     pi.fit.args = pi.fit.args,
                     mu.fit.fun = mu.fit.fun,
                     mu.fit.args = mu.fit.args)
    pi.coef <- coef(res$fit.pi, s = "lambda.min")
    mu.coef <- coef(res$fit.mu, s = "lambda.min")

    pi.vi <- as.numeric(pi.coef != 0)[-1]
    mu.vi <- as.numeric(mu.coef != 0)[-1]
    
    pi.df <- sum(pi.vi)
    mu.df <- sum(mu.vi)
    df <- pi.df + mu.df

    res$fit.pi <- NULL
    res$fit.mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)    
    return(res)
}
