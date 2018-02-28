################################################################
## A class for distributions of p-values:
## log f(p; mu) = g(p)(eta(mu) - eta(mu*)) - (A(mu) - A(mu*)).
## See Section 4.
## Used for deriving the E-step and estimating local FDR.
##
## Attributes:
##    g: transformation of p-value.
##    g.inv: inverse function of g.
##    eta: name of link function. Support all link functions in glm. See ?family for details.
##    mu.star: the parameter value for U([0, 1]).
##    A: log-partition function.
##    name (optional): name of the instance.
##    family (optional): the exponential family used for GLM/GAM. The family for which eta is the canonical link is recommended.
################################################################

exp.family <- function(g, g.inv, eta, mu.star, A,
                       name = NULL, family = NULL){
    f <- function(p, mu){
        ifelse(mu == mu.star, rep(1, length(p)),
               exp(g(p) * (eta(mu) - eta(mu.star)) -
                       (A(mu) - A(mu.star))))
    }
    result <- list(f = f, g = g, g.inv = g.inv, eta = eta,
                   mu.star = mu.star, A = A,
                   name = name, family = family)
    class(result) <- "exp_family"
    return(result)
}

beta.family <- function(){
    g <- function(x){
        tmp <- -log(x)
        pmax(pmin(tmp, -log(10^-15)), -log(1-10^-15))
    }

    g.inv <- function(x){
        tmp <- exp(-x)
        pmax(pmin(tmp, 1 - 10^-15), 10^-15)
    }
    eta <- function(mu){-1/mu}
    mu.star <- 1
    A <- log
    name <- "beta"
    family <- Gamma()
    
    exp.family(g, g.inv, eta, mu.star, A, name, family)
}

inv.gaussian.family <- function(){
    g <- function(x){
        tmp <- qnorm(1 - x)
        pmax(pmin(tmp, qnorm(1 - 10^-15)), qnorm(10^-15))
    }
    g.inv <- function(x){
        tmp <- 1 - pnorm(x)
        pmax(pmin(tmp, 1 - 10^-15), 10^-15)        
    }
    eta <- function(mu){mu}
    mu.star <- 0
    A <- function(mu){mu^2/2}
    name <- "inv_gaussian"
    family <- gaussian()

    exp.family(g, g.inv, eta, mu.star, A, name, family)
}
