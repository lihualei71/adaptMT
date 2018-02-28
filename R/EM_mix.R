################################################################
## EM algorithm to fit a mixture model.
## Combine init_EM_mix.R, Estep.R and Mstep.R together.
##
## Required Input:
##    x: covariate
##    pvals: p-values.
##    plow: lower threshold, i.e. s(x) in AdaPT.
##    phigh: upper threshold, i.e. 1 - s(x) in AdaPT.
##    dist: an "exp_family" object; see exp_family.R.
##    algo: an "AdaPT_algo" object; see gen_EM_algo.R.
##    params0: initial pi(x) and mu(x). Use initilization methods in init_EM_mix.R if pix=NULL or mux=NULL.
##    num.steps: maximal number of steps.
##    tol: tolerance for early stopping. Stop the procedure if ||mux.new - mux.old||_{\infty} < tol.
##    verbose: display a progress bar monitoring the process if verbose = TRUE.
## Output:
##    params: fitted parameters including pix and mux.
################################################################

EM.loglik <- function(dist, params, Hhat, phat){
    pix <- params$pix
    mux <- params$mux
    loglik1 <- sum(Hhat * log(pix) + (1 - Hhat) * log(1 - pix))
    loglik2 <- sum(Hhat * log(dist$f(phat, mux)))
    return(loglik1 + loglik2)
}

EM.mix <- function(x, pvals, plow, phigh, dist,
                   algo = NULL,
                   params0 = list(pix = NULL, mux = NULL),
                   num.steps = 10, tol = 1e-4,
                   verbose = FALSE){
    if (is.null(algo)){
        stop("\"algo\" must be specified! Set algo = list(Mstep.fun = , Mstep.args = , init.fun = , init.args = ).")
    }

    Mstep.fun <- algo$Mstep.fun
    Mstep.args <- algo$Mstep.args
    init.fun <- algo$init.fun
    init.args <- algo$init.args
    input.type <- algo$input.type

    if (verbose){
        cat("Model fitting starts!\n")
        cat("Reduce num.steps if it is too time-consuming!\n")
        pb <- txtProgressBar(min = 0, max = num.steps, style = 3, width = 50)
    }

    if (is.null(params0$pix) || is.null(params0$mux)){
        if (is.null(init.args)){
            init.args <- Mstep.args
        }
        init.args <- c(list(x = x, pvals = pvals, plow = plow,
                            phigh = phigh, dist = dist),
                       init.args)
        result.init <- do.call(init.fun, init.args)
        old.pix <- pix <- result.init$pix
        old.mux <- mux <- result.init$mux
    } else {
        old.pix <- pix <- params0$pix
        old.mux <- mux <- params0$mux
    }

    other.info.list <- list()
    
    for (step in 1:num.steps){
        result.Estep <-
            Estep.mix(x = x, pvals = pvals,
                      plow = plow, phigh = phigh,
                      dist = dist, pix = pix, mux = mux)
        cur.Mstep.args <- c(list(x = x, dist = dist,
                                 Hhat = result.Estep$Hhat,
                                 phat = result.Estep$phat,
                                 input.type = input.type),
                            Mstep.args)
        result.Mstep <- do.call(Mstep.fun, cur.Mstep.args)
        pix <- result.Mstep$pix
        mux <- result.Mstep$mux
        other.info.list[[step]] <- result.Mstep$other
        if (max(abs(mux - old.mux)) < tol &&
            max(abs(pix - old.pix)) < tol){
            break
        }
        old.pix <- pix
        old.mux <- mux
        if (verbose){
            setTxtProgressBar(pb, step)
        }        
    }
    if (verbose){    
        cat("\n")
    }
    params <- list(pix = pix, mux = mux)
    loglik <- EM.loglik(dist, params,
                        result.Estep$Hhat, result.Estep$phat)

    df.vec <- sapply(other.info.list, function(info){info$df})
    df <- mean(df.vec, na.rm = TRUE)

    pi.vi.mat <- sapply(other.info.list, function(info){info$pi.vi})
    if (class(pi.vi.mat)[1] == "matrix"){
        pi.vi <- apply(pi.vi.mat, 1, mean, na.rm = TRUE)
    } else {
        pi.vi <- mean(pi.vi.mat, na.rm = TRUE)
    }

    mu.vi.mat <- sapply(other.info.list, function(info){info$mu.vi})
    if (class(mu.vi.mat)[1] == "matrix"){
        mu.vi <- apply(mu.vi.mat, 1, mean, na.rm = TRUE)
    } else {
        mu.vi <- mean(mu.vi.mat, na.rm = TRUE)
    }
   
    other <- list(df = df, pi.vi = pi.vi, mu.vi = mu.vi)
    
    return(list(params = params, loglik = loglik, other = other))
}

EM.mix.glm <- function(x, pvals, plow, phigh, dist,
                       pi.formula, mu.formula,
                       pi.alter.formulas = NULL,
                       mu.alter.formulas = NULL,
                       params0 = list(pix = NULL, mux = NULL),
                       num.steps = 10, tol = 1e-4,
                       ...){
    algo <- gen.AdaPT.glm(pi.formula = pi.formula,
                          mu.formula = mu.formula,
                          pi.alter.formulas = pi.alter.formulas,
                          mu.alter.formulas = mu.alter.formulas,
                          ...)

    EM.mix(x = x, pvals = pvals, plow = plow, phigh = phigh,
           dist = dist,
           algo = algo, params0 = params0,
           num.steps = num.steps, tol = tol)
}

EM.mix.gam <- function(x, pvals, plow, phigh, dist,
                       pi.formula, mu.formula,
                       pi.alter.formulas = NULL,
                       mu.alter.formulas = NULL,
                       params0 = list(pix = NULL, mux = NULL),
                       num.steps = 10, tol = 1e-4,
                       ...){
    algo <- gen.AdaPT.gam(pi.formula = pi.formula,
                          mu.formula = mu.formula,
                          pi.alter.formulas = pi.alter.formulas,
                          mu.alter.formulas = mu.alter.formulas,
                          ...)

    EM.mix(x = x, pvals = pvals, plow = plow, phigh = phigh,
           dist = dist,
           algo = algo, params0 = params0,
           num.steps = num.steps, tol = tol)
}

EM.mix.glmnet <- function(x, pvals, plow, phigh, dist,
                          params0 = list(pix = NULL, mux = NULL),
                          num.steps = 10, tol = 1e-4,
                          ...){
    algo <- gen.AdaPT.glmnet(...)

    EM.mix(x = x, pvals = pvals, plow = plow, phigh = phigh,
           dist = dist,
           algo = algo, params0 = params0,
           num.steps = num.steps, tol = tol)
}
