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

EM_mix_ms <- function(x, pvals, s, dist, cr = "BIC",
                      cand_methods = NULL,
                      cand_algos = NULL,
                      cand_piargs = NULL,
                      cand_muargs = NULL,
                      params0 = list(pix = NULL, mux = NULL),
                      niter = 20, tol = 1e-4,
                      verbose = TRUE){
    if (is.null(cand_methods) && is.null(cand_algos)){
        stop("Either \"cand_methods\" or \"cand_algos\" must be specified.")
    }

    if (is.null(cand_methods)){
        cand_methods <- rep("custom", length(cand_algos))
    }

    lens <- c(length(cand_methods), length(cand_piargs),
              length(cand_muargs))
    if (max(lens) > min(lens)){
        stop("Number of methods, piargs and muargs should be equal.")
    }
    
    n <- length(pvals)
    info_cr_val <- -Inf
    params <- NULL
    best_model <- NULL
    best_vi <- list(pi = NA, mu = NA)
    best_df <- list(pi = NA, mu = NA)
    
    m <- length(cand_methods)
    if (verbose){
        cat("Model selection starts!\n")
        cat("Shrink the set of candidate models if it is too time-consuming.")
        cat("\n")
        pb <- txtProgressBar(min = 0, max = m, style = 3, width = 50)
    }
    for (i in 1:m){
        method <- cand_methods[i]
        if (method == "custom"){
            algo <- cand_algos[[i]]
        } else {
            algo <- NULL
        }
        piargs <- cand_piargs[[i]]
        muargs <- cand_muargs[[i]]
        mod <- try(
            EM_mix(x, pvals, s, dist, method, algo,
                   piargs, muargs,
                   params0, niter, tol),
            silent = TRUE
            )
        if (class(mod)[1] == "try-error"){
            warning(paste0("Model ", i, " fails."))
            next
        }
        
        loglik <- mod$loglik
        df <- mod$info$pi_df + mod$info$mu_df
        val <- info_cr(loglik, cr, df, n)
        if (val > info_cr_val){
            params <- mod$params
            info_cr_val <- val
            best_model <- list(method = method,
                               algo = algo,
                               piargs = piargs,
                               muargs = muargs)
            best_vi$pi <- mod$info$pi_vi
            best_vi$mu <- mod$info$mu_vi
            best_df$pi <- mod$info$pi_df
            best_df$mu <- mod$info$mu_df
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
    info <- list(pi_df = best_df$pi,
                 mu_df = best_df$mu,
                 pi_vi = best_vi$pi,
                 mu_vi = best_vi$mu)
    return(list(best = best_model,
                params = params,
                info = info))
}
