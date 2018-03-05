#===============================================================
# adapt_model class
#===============================================================

gen_adapt_model_root <- function(pifun, mufun,
                                 pifun_init, mufun_init,
                                 piargs = list(),
                                 muargs = list(),
                                 piargs_init = list(),
                                 muargs_init = list()){
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
        list(algo = algo, args = args),
        class = "adapt_model"
        )
    return(model)
}

gen_adapt_model_glm <- function(pi_formulas, mu_formulas,
                                piargs = NULL, muargs = NULL){
    pifun <- function(formula, data, ...){
        mod <- adapt_glm(formula, data, family = )
    }
}
