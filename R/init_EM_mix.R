#===============================================================
# Initialization of EM
#===============================================================

init_mix_pi <- function(x, pvals, s, fun, args){
    J <- ifelse(
        pvals < s | pvals > 1 - s, 1,
        2 * s / (2 * s - 1)
        )
    J[s == 0.5] <- 0
    
    args <- complete_args(x, J, fun, args)

    if (is.null(args)) {
        stop("pifun_init has irregular input types. Replace another function or writing a wrapper of pifun with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)")
    }

    fit <- do.call(fun, args)
    if (!"fitv" %in% names(fit)){
        stop("pifun_init does not output fitv. Replace another function or change the name for fitted value to fitv")
    }

    pix <- as.numeric(fit$fitv)
    if (any(is.na(pix))){
        stop("Initialization of pix has NAs")
    }
    pix <- pminmax(pix, 0, 1)
    
    return(
        list(pix = pix,
             fit_pi = fit)
        )
}

init_mix_mu <- function(x, pvals, s, dist, fun, args){
    imputed_p <- ifelse(
        pvals < s | pvals > 1 - s,
        pmin(pvals, 1 - pvals),
        pvals
        )
    imputed_p <- pminmax(imputed_p, 1e-15, 1-1e-15)
    
    args <- complete_args(x, J, fun, args)

    if (is.null(args)) {
        stop("mufun_init has irregular input types. Replace another function or writing a wrapper of mufun with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)")
    }

    fit <- do.call(fun, args)
    if (!"fitv" %in% names(fit)){
        stop("mufun_init does not output fitv. Replace another function or change the name for fitted value to fitv")
    }

    mux <- as.numeric(fit$fitv)
    if (any(is.na(mux))){
        stop("Initialization of mux has NAs")
    }
    if (dist$family$family == "Gamma"){
        mux <- pmax(mux, 1)
    } else if (dist$family$family == "gaussian"){
        mux <- pmax(mux, 0)
    }
    
    return(
        list(mux = mux,
             fit_mu = fit)
        )
}

init_mix_root <- function(x, pvals, s, dist,
                          pifun, mufun,
                          piargs = NULL, muargs = NULL){
    pi_res <- init_mix_pi(x, pvals, s, pifun, piargs)
    mu_res <- init_mix_mu(x, pvals, s, dist, mufun, muargs)
    res <- c(pi_res, mu_res)
    pi_info <- modinfo(res$fit_pi)
    mu_info <- modinfo(res$fit_mu)
    res$info <- c(pi_info, mu_info)
    res$fit_pi <- res$fit_mu <- NULL
    return(res)
}

init_mix_glm <- function(x, pvals, s, dist,
                         pi_formula, mu_formula,
                         piargs = NULL, muargs = NULL){
    pifun <- function(formula, data, ...){
        safe.glm(formula, data, 
                 family = gaussian(), ...)
    }
    piargs <- c(list(formula = pi_formula), piargs)

    mufun <- function(formula, data, weights, ...){
        safe.glm(formula, data, weights = weights,
                 family = dist$family, ...)
    }
    muargs <- c(list(formula = mu_formula), muargs)

    res <- init_mix_root(x, pvals, s, dist,
                         pifun, mufun,
                         piargs, muargs)
    return(res)    
}

init_mix_gam <- function(x, pvals, s, dist,
                         pi_formula, mu_formula,
                         piargs = NULL, muargs = NULL){
    pifun <- function(formula, data, ...){
        safe.gam(formula, data, 
                 family = gaussian(), ...)
    }
    piargs <- c(list(formula = pi_formula), piargs)

    mufun <- function(formula, data, weights, ...){
        safe.gam(formula, data, weights = weights,
                 family = dist$family, ...)
    }
    muargs <- c(list(formula = mu_formula), muargs)

    res <- init_mix_root(x, pvals, s, dist,
                         pifun, mufun,
                         piargs, muargs)
    return(res)
}

init_mix_glmnet <- function(x, pvals, s, dist,
                            piargs, muargs){
    pifun <- function(x, y, ...){
        safe.glmnet(x, y, 
                    family = "gaussian", ...)
    }
    mufun <- function(x, y, weights, ...){
        safe.glmnet(x, y, weights = weights,
                    family = dist$family, ...)
    }

    res <- init_mix_root(x, pvals, s, dist,
                         pifun, mufun,
                         piargs, muargs)
    return(res)
}

init_mix <- function(x, pvals, s, dist,
                     method = c("glm", "gam", "glmnet", "custom"),
                     pifun = NULL, mufun = NULL,
                     piargs = NULL, muargs = NULL,
                     ...){
    method <- method[1]
    func <- switch(
        method,
        glm = init_mix_glm,
        gam = init_mix_gam,
        glmnet = init_mix_glmnet,
        custom = NA)

    if (!is.na(func)){
        func(x, pvals, s, dist,
             piargs = piargs, muargs = muargs,
             ...)
    } else {
        init_mix_root(x, pvals, s, dist,
                      pifun, mufun,
                      piargs, muargs)
    }
}
