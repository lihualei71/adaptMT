################################################################
## Model selection based on partially masked data using information criteria.
##
## Required Input:
##    x: covariate
##    p: p-values.
##    plow: lower threshold, i.e. s(x) in AdaPT.
##    phigh: upper threshold, i.e. 1 - s(x) in AdaPT.
##    dist: distribution family for p-values in "exp_family" class.
##    cand.args: a list of arguments passed into EM.mix. Each represent a candidate model.
##    cr: AIC or BIC.
##    ...: other arguments passed into EM.mix().
## Output:
##    model: selected model. an element from cv.args.
################################################################

info.cr <- function(loglik, cr, df, n){
    switch(cr,
           "AIC" = 2 * df - 2 * loglik,
           "AICC" = 2 * df * n / (n - df - 1),
           "BIC" = log(n) * df - 2 * loglik,
           "HIC" = 2 * log(log(n)) * df - 2 * loglik)
}

EM.mix.ms <- function(x, pvals, plow, phigh, dist,
                      cand.algos = NULL,
                      cr = "BIC",
                      params0 = list(pix = NULL, mux = NULL),
                      num.steps = 20, tol = 1e-4,
                      verbose = TRUE){
    if (is.null(cand.algos)){
        stop("cand.algos must be specified! Set cand.algos as a list of \"AdaPT_algo\" objects.")
    }
    
    n <- length(pvals)
    info.cr.val <- -Inf
    params <- NULL
    best.model <- NULL
    best.vi <- list(pi = NA, mu = NA)
    best.df <- NA
    
    m <- length(cand.algos)
    if (verbose){
        cat("Model selection starts!\n")
        cat("Shrink the set of candidate models if it is too time-consuming!\n")
        pb <- txtProgressBar(min = 0, max = m, style = 3, width = 50)
    }
    for (i in 1:m){
        algo <- cand.algos[[i]]
        mod <- try(EM.mix(x = x, pvals = pvals,
                          plow = plow, phigh = phigh, dist = dist,
                          algo = algo, params0 = params0,
                          num.steps = num.steps, tol = tol))
        if (class(mod)[1] == "try-error"){
            next
        }
        
        loglik <- mod$loglik
        df <- mod$other$df
        val <- info.cr(loglik, cr, df, n)
        if (val > info.cr.val){
            params <- mod$params
            info.cr.val <- val
            best.model <- algo
            best.vi$pi <- mod$pi.vi
            best.vi$mu <- mod$mu.vi
            best.df <- df
        }

        if (verbose){
            setTxtProgressBar(pb, i)
        }
    }
    if (verbose){    
        cat("\n")
    }
    if (info.cr.val == -Inf){
        stop("All Models Fail!")
    }
    other <- list(df = best.df,
                  pi.vi = best.vi$pi,
                  mu.vi = best.vi$mu)
    return(list(best = best.model,
                params = params,
                other = other))
}

