#---------------------------------------------------------------
# Model selection based on partially masked data
#---------------------------------------------------------------

info_cr <- function(loglik, cr, df, n){
    switch(cr,
           "AIC" = 2 * df - 2 * loglik,
           "AICC" = 2 * df * n / (n - df - 1),
           "BIC" = log(n) * df - 2 * loglik,
           "HIC" = 2 * log(log(n)) * df - 2 * loglik)
}

EM_mix_ms <- function(x, pvals, s, dist, models,
                      cr = "BIC",
                      params0 = list(pix = NULL, mux = NULL),
                      niter = 20, tol = 1e-4,
                      verbose = TRUE,
                      type = "unweighted"){
    n <- length(pvals)
    info_cr_val <- -Inf
    
    m <- length(models)
    if (verbose){
        cat("Model selection starts!\n")
        cat("Shrink the set of candidate models if it is too time-consuming.")
        cat("\n")
        pb <- txtProgressBar(min = 0, max = m, style = 3, width = 50)
    }
    for (i in 1:m){
        model <- complete_model(models[[i]], dist)
        fit <- try(
            EM_mix(x, pvals, s, dist, model, params0, niter, tol,
                   type = type),
            silent = TRUE
            )
        if (class(fit)[1] == "try-error"){
            warning(paste0("Model ", i, " fails."))
            next
        }
        
        loglik <- fit$loglik
        df <- fit$info$pi$df + fit$info$mu$df
        val <- info_cr(loglik, cr, df, n)
        if (val > info_cr_val){
            params <- fit$params
            info_cr_val <- val
            best_model <- models[[i]]
            best_model_info <- fit$info
        }

        if (verbose){
            setTxtProgressBar(pb, i)
        }
    }
    if (verbose){    
        cat("\n")
    }
    if (info_cr_val == -Inf){
        stop("All models fail.")
    }
    return(list(model = best_model,
                params = params,
                info = best_model_info))
}
