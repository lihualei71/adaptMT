#===============================================================
# adapt_model class
#===============================================================

gen_adapt_model <- function(pifun, mufun,
                            pifun_init, mufun_init,
                            piargs = list(),
                            muargs = list(),
                            piargs_init = list(),
                            muargs_init = list(),
                            name = ""){
    if (!is.function(pifun)){
        stop("\"pifun\" must be a function.")
    }
    if (!is.function(mufun)){
        stop("\"mufun\" must be a function.")
    }
    if (!is.function(pifun_init)){
        stop("\"pifun_init\" must be a function.")
    }
    if (!is.function(mufun_init)){
        stop("\"mufun_init\" must be a function.")
    }

    algo <- list(pifun = pifun, mufun = mufun,
                 pifun_init = pifun_init, mufun_init = mufun_init)
    args <- list(piargs = piargs, muargs = muargs,
                 piargs_init = piargs_init, muargs_init = muargs_init)
    model <- structure(
        list(name = name, algo = algo, args = args),
        class = "adapt_model"
        )
    return(model)
}

gen_adapt_model_glm <- function(dist, 
                                piargs = list(),
                                muargs = list()){
    pifun <- function(formula, data, ...){
        adapt_glm(formula, data,
                  family = binomial(), ...)
    }
    
    pifun_init <- function(formula, data, ...){
        adapt_glm(formula, data,
                  family = gaussian(), ...)
    }
    
    mufun <- mufun_init <- function(formula, data, weights, ...){
        adapt_glm(formula, data, weights = weights,
                  family = dist$family, ...)
    }

    piargs_init <- piargs
    muargs_init <- muargs

    gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
                    piargs, muargs, piargs_init, muargs_init)
}

gen_adapt_model_gam <- function(dist, 
                                piargs = list(),
                                muargs = list()){
    pifun <- function(formula, data, ...){
        adapt_gam(formula, data,
                  family = binomial(), ...)
    }
    
    pifun_init <- function(formula, data, ...){
        adapt_gam(formula, data,
                  family = gaussian(), ...)
    }
    
    mufun <- mufun_init <- function(formula, data, weights, ...){
        adapt_gam(formula, data, weights = weights,
                  family = dist$family, ...)
    }

    piargs_init <- piargs
    muargs_init <- muargs

    gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
                    piargs, muargs, piargs_init, muargs_init)
}

gen_adapt_model_glmnet <- function(dist, 
                                   piargs = list(),
                                   muargs = list()){
    pifun <- function(x, y, ...){
        adapt_glmnet(x, y,
                     family = "binomial", ...)
    }
    
    pifun_init <- function(x, y, ...){
        adapt_glmnet(x, y,
                     family = "gaussian", ...)
    }
    
    mufun <- mufun_init <- function(x, y, weights, ...){
        adapt_glmnet(x, y, weights = weights,
                     family = dist$family, ...)
    }

    piargs_init <- piargs
    muargs_init <- muargs

    gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
                    piargs, muargs, piargs_init, muargs_init)
}
