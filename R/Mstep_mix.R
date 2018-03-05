#===============================================================
# Compute M-step for the mixture model.
#===============================================================

Mstep_mix_pi <- function(x, Hhat, fun, args){
    args <- complete_args(x, Hhat, fun, args)

    if (is.null(args)) {
        stop("pifun has irregular input types. Replace another function or writing a wrapper of pifun with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)")
    }

    fit <- do.call(fun, args)
    if (is.null(fit$fitv)){
        stop("pifun does not output fitv. Replace another function or change the name for fitted value to fitv")
    }

    pix <- as.numeric(fit$fitv)
    if (any(is.na(pix))){
        stop("pix in M-step has NAs")
    }
    pix <- pminmax(pix, 0, 1)
    
    return(
        list(pix = pix,
             fit_pi = fit$mod)
        )
}

Mstep_mix_mu <- function(x, Hhat, phat, dist, fun, args){
    yhat <- dist$g(phat)
    Hhat <- pmax(Hhat, 1e-5)
    if (!"weights" %in% formalArgs(fun)){
        stop("mufun does not have input weights")
    }

    args <- complete_args(x, yhat, fun, args, Hhat)

    if (is.null(args)){
        stop("mufun has irregular input types. Replace another function or writing a wrapper of mufun with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)")
    }

    fit <- do.call(fun, args)
    if (is.null(fit$fitv)){
        stop("mufun does not output fitv. Replace another function or change the name for fitted value to fitv")
    }

    mux <- as.numeric(fit$fitv)
    if (any(is.na(mux))){
        stop("mux in M-step has NAs")
    }
    if (dist$family$family == "Gamma"){
        mux <- pmax(mux, 1)
    } else if (dist$family$family == "gaussian"){
        mux <- pmax(mux, 0)
    }
    
    return(
        list(mux = mux,
             fit_mu = fit$mod)
        )
}

Mstep_mix_root <- function(x, Hhat, phat, dist,
                           pifun, mufun, 
                           piargs = NULL, muargs = NULL){
    pi_res <- Mstep_mix_pi(x = x, Hhat = Hhat,
                           fun = pifun, args = piargs)
    mu_res <- Mstep_mix_mu(x = x, Hhat = Hhat, phat = phat,
                           dist = dist,
                           fun = mufun, args = muargs)
    res <- c(pi_res, mu_res)
    pi_info <- modinfo(res$fit_pi)
    mu_info <- modinfo(res$fit_mu)
    res$info <- list(pi_df = pi_info$df, pi_vi = pi_info$vi,
                     mu_df = mu_info$df, mu_vi = mu_info$vi)
    res$fit_pi <- res$fit_mu <- NULL
    return(res)
}

Mstep_mix_glm <- function(x, Hhat, phat, dist,
                          piargs = NULL, muargs = NULL){
    pifun <- function(formula, data, ...){
        adapt_glm(formula, data,
                 family = binomial(), ...)
    }
    if (is.null(piargs$formula)){
        stop("argument \"formula\" is missing. Please specify it in \"piargs\"")
    }

    mufun <- function(formula, data, weights, ...){
        adapt_glm(formula, data, weights = weights,
                 family = dist$family, ...)
    }
    if (is.null(muargs$formula)){
        stop("argument \"formula\" is missing. Please specify it in \"muargs\"")
    }

    res <- Mstep_mix_root(x, Hhat, phat, dist,
                          pifun, mufun,
                          piargs, muargs)
    return(res)
}

Mstep_mix_gam <- function(x, Hhat, phat, dist,
                          piargs = NULL, muargs = NULL){
    pifun <- function(formula, data, ...){
        adapt_gam(formula, data,
                 family = binomial(), ...)
    }
    if (is.null(piargs$formula)){
        stop("argument \"formula\" is missing. Please specify it in \"piargs\"")
    }
    
    mufun <- function(formula, data, weights, ...){
        adapt_gam(formula, data, weights = weights,
                 family = dist$family, ...)
    }
    if (is.null(muargs$formula)){
        stop("argument \"formula\" is missing. Please specify it in \"muargs\"")
    }

    res <- Mstep_mix_root(x, Hhat, phat, dist,
                          pifun, mufun,
                          piargs, muargs)
    return(res)
}

Mstep_mix_glmnet <- function(x, Hhat, phat, dist,
                             piargs = NULL, muargs = NULL){
    pifun <- function(x, y, ...){
        adapt_glmnet(x, y,
                    family = "binomial", ...)
    }

    mufun <- function(x, y, weights, ...){
        adapt_glmnet(x, y, weights = weights,
                    family = dist$family, ...)
    }

    res <- Mstep_mix_root(x, Hhat, phat, dist,
                          pifun, mufun,
                          piargs, muargs)
    return(res)
}

Mstep_mix <- function(x, Hhat, phat, dist,
                      method = c("glm", "gam", "glmnet", "custom"),
                      pifun = NULL, mufun = NULL,
                      piargs = NULL, muargs = NULL){
    method <- method[1]
    func <- switch(
        method,
        glm = Mstep_mix_glm,
        gam = Mstep_mix_gam,
        glmnet = Mstep_mix_glmnet,
        custom = NA)

    if (method != "custom"){
        func(x, Hhat, phat, dist,
             piargs = piargs, muargs = muargs)
    } else {
        Mstep_mix_root(x, Hhat, phat, dist,
                       pifun, mufun,
                       piargs, muargs)
    }
}
