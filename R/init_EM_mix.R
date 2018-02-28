################################################################
## Parameter initialization for the mixture model.
## 
## Required Input:
##    x: covariate.
##    pvals: p-values.
##    plow: lower threshold, i.e. s0(x) in AdaPT.
##    phigh: upper threshold, i.e. 1 - s0(x) in AdaPT.
##    dist: distribution family for p-values in "exp_family" class.
##    pi.fit.fun: a function to fit initial pix using \tilde{J_i}.
##    pi.fit.args: other arguments passed into pi.fit.fun.
##    mu.fit.fun: a model to fit initial mux using (x,min(p,1-p)).
##    mu.fit.args: other arguments passed into mu.fit.fun.
## Output:
##    pix, mux: initial guess of pi(x) and mu(x).
################################################################

init.mix <- function(x, pvals, plow, phigh, dist,
                     input.type = c("formula", "xy"),
                     pi.fit.fun = NULL,
                     pi.fit.args = NULL,
                     mu.fit.fun = NULL,
                     mu.fit.args = NULL){
    input.type <- input.type[1]

    J <- ifelse(pvals < plow | pvals > phigh, 1,
                (plow + 1 - phigh) / (plow - phigh))
    imputed.p <- ifelse(pvals < plow | pvals > phigh,
                        pmin(pvals, 1 - pvals), pvals)
    imputed.p <- pmin(pmax(imputed.p, 1e-15), 1-1e-15)

    Mstep.mix(x = x, Hhat = J, phat = imputed.p, dist = dist,
              input.type = input.type,
              pi.fit.fun = pi.fit.fun,
              pi.fit.args = pi.fit.args,
              mu.fit.fun = mu.fit.fun,
              mu.fit.args = mu.fit.args)
}

init.mix.glm <- function(x, pvals, plow, phigh, dist,
                         pi.formula, mu.formula,
                         input.type = "formula",
                         ...){
    if (!is.null(input.type) && (input.type != "formula")){
        input.type <- "formula"
        warning("Warning: input.type for GLM can only be formula!")
    }
    extra.args <- list(...)
    pi.fit.fun <- safe.gaussian.glm
    pi.fit.args <- c(list(formula = pi.formula), extra.args)
    mu.fit.fun <- function(formula, data, ...){
        safe.glm(formula, data, 
                 family = dist$family, ...)
    }
    mu.fit.args <- c(list(formula = mu.formula), extra.args)
    res <- init.mix(x = x, pvals = pvals,
                    plow = plow, phigh = phigh,
                    dist = dist,
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

init.mix.gam <- function(x, pvals, plow, phigh, dist,
                         pi.formula, mu.formula,
                         input.type = "formula",
                         ...){
    if (!is.null(input.type) && (input.type != "formula")){
        input.type <- "formula"
        warning("Warning: input.type for GLM can only be formula!")
    }
    extra.args <- list(...)    
    pi.fit.fun <- safe.gaussian.gam
    pi.fit.args <- c(list(formula = pi.formula), extra.args)
    mu.fit.fun <- function(formula, data, ...){
        safe.gam(formula, data, 
                 family = dist$family, ...)
    }
    mu.fit.args <- c(list(formula = mu.formula), extra.args)
    res <- init.mix(x = x, pvals = pvals,
                    plow = plow, phigh = phigh,
                    dist = dist,
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

init.mix.glmnet <- function(x, pvals, plow, phigh, dist,
                            input.type = "xy",
                            ...){
    if (!is.null(input.type) && (input.type != "xy")){
        input.type <- "xy"
        warning("Warning: input.type for glmnet can only be xy!")
    }
    extra.args <- list(...)
    pi.fit.fun <- safe.gaussian.glmnet
    pi.fit.args <- extra.args
    mu.fit.fun <- function(x, y, ...){
        safe.glmnet(x, y, 
                    family = dist$family, ...)
    }
    mu.fit.args <- extra.args
    res <- init.mix(x = x, pvals = pvals,
                    plow = plow, phigh = phigh,
                    dist = dist,
                    input.type = input.type,
                    pi.fit.fun = pi.fit.fun,
                    pi.fit.args = pi.fit.args,
                    mu.fit.fun = mu.fit.fun,
                    mu.fit.args = mu.fit.args)
    pi.coef <- coef(res$fit.pi, s = "lambda.min")
    mu.coef <- coef(res$fit.mu, s = "lambda.min")

    pi.vi <- as.numeric(pi.coef != 0)
    mu.vi <- as.numeric(mu.coef != 0)
    
    pi.df <- sum(pi.vi)
    mu.df <- sum(mu.vi)
    df <- pi.df + mu.df

    res$fit.pi <- NULL
    res$fit.mu <- NULL
    res$other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)
    return(res)    
}
