##############################################################################
## AdaPT with varying-coefficient two-group model
## Lihua Lei & William Fithian,
##   "AdaPT: An interactive procedure for multiple testing with side information"
## Available from http://arxiv.org/abs/1609.06035
##############################################################################

find.snew.mix <- function(x, pvals, s, dist,
                          params, delta){
    mask <- (pvals <= s) | (pvals >= 1 - s)
    R <- sum(mask)
    mask.pvals <- pmin(pvals, 1 - pvals)
    mask.lfdr <- lfdr.mix(mask.pvals, dist, params)
    sorted.lfdr.masked <- sort(mask.lfdr[mask], decreasing = TRUE)
    num.reveal <- ceiling(R * delta)
    lfdr.lev <- (sorted.lfdr.masked[num.reveal] +
                     sorted.lfdr.masked[num.reveal + 1]) / 2
    lfdr.lev <- as.numeric(lfdr.lev)
    k <- 1
    
    while (TRUE){
        s.new <- compute.threshold.mix(dist, params, lfdr.lev)
        s.new <- pmin(s, s.new)
        R.new <- sum((pvals <= s.new) | (pvals >= 1-s.new))
        if (R.new <= R - num.reveal){
            return(s.new)
        } else {
            lfdr.lev <- lfdr.lev * (1 - 1e-4 * k)
            k <- k + 1
        }
    }
    return(s.new)
}

AdaPT <- function(x, pvals, dist,
                  algo = NULL,
                  EM.args = NULL,
                  ms.fun = EM.mix.ms, cand.algos = NULL,
                  ms.args = NULL,
                  plot.fun = NULL, plot.args = NULL,
                  s0 = rep(0.45, length(pvals)),
                  params0 = list(pix = NULL, mux = NULL),
                  delta = 1/length(pvals), R.min = 1,
                  delta.up = 0.05, delta.ms = 1,
                  alphas = seq(0.01, 1, 0.01),
                  print.verbose = TRUE,
                  up.verbose = FALSE,
                  ms.verbose = TRUE,
                  return.data = TRUE){
    ## Record main information
    n <- length(pvals) # number of hypotheses
    m <- length(alphas) # number of alpha's
    s.return <- matrix(0, n, m) # threshold
    fdp.return <- rep(0, m) # post-hoc fdp estimate for each alpha
    num_rej.return <- rep(0, m) # number of rejections
    params.return <- list() # parameters (including pix and mux)
    lfdr.return <- matrix(0, n, m) # local fdr estimate
    algo.list <- list() # all selected models
    pi.vi.list <- list() # variable importance of pi
    mu.vi.list <- list() # variable importance of mu
    
    ## Initialization
    params <- params0
    R.min <- max(R.min, sum(pvals <= 0) + sum(pvals >= 1))
    fit.args.root <- c(list(x = x, pvals = pvals, dist = dist,
                            verbose = up.verbose),
                       EM.args)
    ## Perform model selection if "algo" is not specified and "cand.algos" is not specified
    ms.flag <- !is.null(cand.algos) && is.null(algo)
    if (ms.flag){
        ms.args.root <- c(list(x = x, pvals = pvals, dist = dist,
                               cand.algos = cand.algos,
                               verbose = ms.verbose),
                          ms.args)
    } else {
        algo.list[[1]] <- algo
    }
    if (!is.null(plot.fun)){
        plot.args.root <- c(list(x = x), plot.args)
    }
    s <- s0
    fdp <- fdp.hat(pvals, s)
    alphaind <- m
    alpha <- alphas[alphaind]
    step <- 0
    period.up <- ceiling(delta.up * n)
    period.ms <- ceiling(delta.ms * n)
    
    while (fdp > alpha | alphaind > 0){
        if (ms.flag && (step %% period.ms == 0)){
            ms.args <- c(list(plow = s, phigh = 1 - s,
                              params0 = params),
                         ms.args.root)
            ms.res <- do.call(ms.fun, ms.args)
            algo <- ms.res$best
            params <- ms.res$params
            algo.list <- append(algo.list, list(algo))
            pi.vi.list <- append(pi.vi.list,
                                 list(ms.res$other$pi.vi))
            mu.vi.list <- append(mu.vi.list,
                                 list(ms.res$other$mu.vi))
        }
        
        if (step %% period.up == 0){
            fit.args <- c(list(plow = s, phigh = 1 - s,
                               params0 = params,
                               algo = algo),
                          fit.args.root)
            mod <- do.call(EM.mix, fit.args)
            params <- mod$params
            lfdr <- lfdr.mix(pvals, dist, params)
            pi.vi.list <- append(pi.vi.list,
                                 list(mod$other$pi.vi))
            mu.vi.list <- append(mu.vi.list,
                                 list(mod$other$mu.vi))
        }

        s.new <- find.snew.mix(x, pvals, s, dist, params, delta)
        s.new <- pmin(s, s.new)

        if (!is.null(plot.fun)){
            plot.args <- c(list(s = s, s.new = s.new,
                                params = params,
                                lfdr = lfdr),
                           plot.args.root)
            do.call(plot.fun, plot.args)
        }

        step <- step + 1        
        s <- s.new
        fdp <- fdp.hat(pvals, s)
        R <- sum(pvals <= s)
        if (R <= R.min) {
            break
        }
        if (fdp <= alpha) {
            if (print.verbose){
                cat(paste0("Step ", step, ": FDP ", fdp, ", Number of Rej. ", R, "\n"))
            }

            remain.alpha <- which(alphas < fdp)
            if (length(remain.alpha) == 0){
                exit.flag <- TRUE
                new.alphaind <- 0
            } else {
                exit.flag <- FALSE
                new.alphaind <- max(remain.alpha)
            }
            for (j in (new.alphaind + 1):alphaind){
                s.return[, j] <- s
                num_rej.return[j] <- R
                fdp.return[j] <- fdp
                params.return[[j]] <- params
                lfdr.return[, j] <- lfdr
            }
            if (exit.flag) break
            alphaind <- new.alphaind
            alpha <- alphas[alphaind]
        }
    }

    pi.vi.mat <- Reduce(cbind, pi.vi.list)
    if (class(pi.vi.mat)[1] == "matrix"){
        pi.vi <- apply(pi.vi.mat, 1, mean, na.rm = TRUE)
    } else {
        pi.vi <- mean(pi.vi.mat, na.rm = TRUE)
    }
    mu.vi.mat <- Reduce(cbind, mu.vi.list)
    if (class(mu.vi.mat)[1] == "matrix"){
        mu.vi <- apply(mu.vi.mat, 1, mean, na.rm = TRUE)
    } else {
        mu.vi <- mean(mu.vi.mat, na.rm = TRUE)
    }
    
    info.list <- list(vi = list(pi = pi.vi, mu = mu.vi),
                      dist = dist, delta = delta,
                      cand.algos = cand.algos,
                      ms.fun = ms.fun,
                      EM = EM.args,
                      ms = ms.args)
    results <- list(num_rej = num_rej.return,
                    fdp = fdp.return,
                    s = s.return,
                    params = params.return,
                    lfdr = lfdr.return,
                    alphas = alphas,
                    algos = algo.list,
                    info = info.list)
    if (return.data){
        results$data <- list(x = x, pvals = pvals)
    }
    class(results) <- c("AdaPT")
    return(invisible(results))
}

AdaPT.glm <- function(x, pvals, dist,
                      pi.formulas, mu.formulas,
                      glm.args = NULL,
                      num.steps.up = 10,
                      num.steps.ms = 30,
                      tol = 1e-4,
                      ...){
    if (class(pi.formulas)[1] == "list"){
        pi.formulas <- unlist(pi.formulas)
    }
    if (class(mu.formulas)[1] == "list"){
        mu.formulas <- unlist(mu.formulas)
    }
    nf <- length(pi.formulas)
    if (nf != length(mu.formulas)){
        stop("the numbers of formulas for pi(x) and mu(x) do not match!")
    }
    if (nf == 1){
        gen.algo.args <- c(list(pi.formula = pi.formulas,
                                mu.formula = mu.formulas),
                           glm.args)
        algo <- do.call(gen.AdaPT.glm, gen.algo.args)
        cand.algos <- NULL
    } else {
        algo <- NULL
        cand.algos <- lapply(1:nf, function(i){
            gen.algo.args <- c(list(pi.formula = pi.formulas[i],
                                    mu.formula = mu.formulas[i]),
                               glm.args)
            do.call(gen.AdaPT.glm, gen.algo.args)
        })
    }
    EM.args <- list(num.steps = num.steps.up, tol = tol)
    ms.args <- list(num.steps = num.steps.ms, tol = tol)
    ms.fun <- EM.mix.ms
    AdaPT(x = x, pvals = pvals, dist = dist,
          algo = algo, EM.args = EM.args,
          ms.fun = ms.fun, cand.algos = cand.algos,
          ms.args = ms.args, ...)
}

AdaPT.gam <- function(x, pvals, dist,
                      pi.formulas, mu.formulas,
                      gam.args = NULL,
                      num.steps.up = 10,
                      num.steps.ms = 30,
                      tol = 1e-4,
                      ...){
    if (class(pi.formulas)[1] == "list"){
        pi.formulas <- unlist(pi.formulas)
    }
    if (class(mu.formulas)[1] == "list"){
        mu.formulas <- unlist(mu.formulas)
    }
    nf <- length(pi.formulas)
    if (nf != length(mu.formulas)){
        stop("Error: the numbers of formulas for pi(x) and mu(x) do not match!")
    }
    if (nf == 1){
        gen.algo.args <- c(list(pi.formula = pi.formulas,
                                mu.formula = mu.formulas),
                           gam.args)
        algo <- do.call(gen.AdaPT.gam, gen.algo.args)
        cand.algos <- NULL
    } else {
        algo <- NULL
        cand.algos <- lapply(1:m, function(i){
            gen.algo.args <- c(list(pi.formula = pi.formulas[i],
                                    mu.formula = mu.formulas[i]),
                               gam.args)
            do.call(gen.AdaPT.gam, gen.algo.args)
        })
    }
    EM.args <- list(num.steps = num.steps.up, tol = tol)
    ms.args <- list(num.steps = num.steps.ms, tol = tol)
    ms.fun <- EM.mix.ms
    AdaPT(x = x, pvals = pvals, dist = dist,
          algo = algo, EM.args = EM.args,
          ms.fun = ms.fun, cand.algos = cand.algos,
          ms.args = ms.args, ...)
}

AdaPT.glmnet <- function(x, pvals, dist,
                         glmnet.args = NULL,
                         num.steps = 10,
                         tol = 1e-4,
                         ...){
    algo <- do.call(gen.AdaPT.glmnet, glmnet.args)
    cand.algos <- NULL
    EM.args <- list(num.steps = num.steps, tol = tol)
    ms.fun <- EM.mix.ms
    AdaPT(x = x, pvals = pvals, dist = dist,
          algo = algo, EM.args = EM.args,
          ms.fun = ms.fun, cand.algos = cand.algos, ...)
}

