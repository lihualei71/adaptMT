################################################################
## Functions to generate common "AdaPT_algo" objects
## An "AdaPT_algo" object is a list containing the model fitting in## formation with the following components:
##    Mstep.fun: function to perform M-step; see Mstep_mix.R;
##    Mstep.args: arguments of Mstep.fun;
##    init.fun: function for initialization; see init_EM_mix.R;
##    init.args: arguments of init.fun;
##    input.type: either "formula" or "xy".
##
## Required Input:
##    name: the name of the algorithm. Supporting "glm", "gam" or "glmnet" in current version;
##    Mstep.args: arguments of Mstep.fun;
##    init.args: arguments of init.fun. Default to be the same as Mstep.args.
## 
## Output:
##    algo: an "AdaPT_algo" object.
################################################################

gen.AdaPT.glm <- function(pi.formula, mu.formula,
                          ...){
    extra.args <- list(...)
    Mstep.args <- init.args <-
        c(list(pi.formula = pi.formula,
               mu.formula = mu.formula),
          extra.args)
    obj <- list(Mstep.fun = Mstep.mix.glm,
                Mstep.args = Mstep.args,
                init.fun = init.mix.glm,
                init.args = init.args,
                input.type = "formula")
    class(obj) <- "AdaPT_algo"
    return(obj)
}

gen.AdaPT.gam <- function(pi.formula, mu.formula,
                          ...){
    extra.args <- list(...)
    Mstep.args <- init.args <-
        c(list(pi.formula = pi.formula,
               mu.formula = mu.formula),
          extra.args)
    obj <- list(Mstep.fun = Mstep.mix.gam,
                Mstep.args = Mstep.args,
                init.fun = init.mix.gam,
                init.args = init.args,
                input.type = "formula")
    class(obj) <- "AdaPT_algo"
    return(obj)
}

gen.AdaPT.glmnet <- function(...){
    Mstep.args <- init.args <- list(...)
    obj <- list(Mstep.fun = Mstep.mix.glmnet,
                Mstep.args = Mstep.args,
                init.fun = init.mix.glmnet,
                init.args = init.args,
                input.type = "xy")
    class(obj) <- "AdaPT_algo"
    return(obj)
}
