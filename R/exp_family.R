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
##    mustar: the parameter value for U([0, 1]).
##    A: log-partition function.
##    name (optional): name of the instance.
##    family (optional): the exponential family used for GLM/GAM. The family for which eta is the canonical link is recommended.
################################################################

#' An S3generic to define an exponential family
#'
#' An object containing all required information in an exponential
#' family to perform the E-step. The exponential function is encoded by
#' \deqn{h(p; \eta) = \exp\{(\eta(\mu) - \eta(\mu^{*})) g(p) - (A(\mu) - A(\mu^{*}))\}}{h(p; \eta) = exp{(\eta(\mu) - \eta(\mu*)) g(p) - (A(\mu) - A(\mu*))}}
#' where \eqn{g(p)} is an arbitrary transformation, \eqn{\mu} is the
#' \emph{mean parameter}, \eqn{\eta} is the natural parameter,
#' and \eqn{A(\mu)} is the partition function. The extra redundant
#' parameter \eqn{\mu^{*}}{\mu*} is to guarantee that \eqn{U([0, 1])}
#' belongs to the class.
#' 
#' @param g function. An transformation of p-values
#' @param ginv function. The inverse function of \code{g}
#' @param eta function. The natural parameter as a function of the mean parameter \code{mu}
#' @param mustar scalar. The mean parameter that gives \eqn{U([0, 1])}
#' @param A function. The partition function
#' @param name character. A name for the family. NULL by default
#' @param family object of S3 class "\code{\link[stats]{gaussian}}" from \code{stats} package. The family used for model fitting in \code{glm}, \code{gam}, \code{glmnet}, etc..
#' @return an object of S3 class "exp_family". This includes all inputs and  \item{h }{function. The density function computed}
#' @family aggregate functions
#' @seealso \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'    sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over
#'    high-dimensional arrays.
#' 
#'
#' @export
#' 
exp_family <- function(g, ginv, eta, mustar, A,
                       name = NULL, family = NULL){
    h <- function(p, mu){
        ifelse(mu == mustar, rep(1, length(p)),
               exp(g(p) * (eta(mu) - eta(mustar)) -
                       (A(mu) - A(mustar))))
    }
    result <- list(h = h, g = g, ginv = ginv, eta = eta,
                   mustar = mustar, A = A,
                   name = name, family = family)
    class(result) <- "exp_family"
    return(result)
}

#' @rdname exp_family 
#' 
#' @export
#' 
beta_family <- function(){
    g <- function(x){
        tmp <- -log(x)
        pmax(pmin(tmp, -log(10^-15)), -log(1-10^-15))
    }

    ginv <- function(x){
        tmp <- exp(-x)
        pmax(pmin(tmp, 1 - 10^-15), 10^-15)
    }
    eta <- function(mu){-1/mu}
    mustar <- 1
    A <- log
    name <- "beta"
    family <- Gamma()

    exp_family(g, ginv, eta, mustar, A, name, family)
}

#' @rdname exp_family 
#' 
#' @export
#' 
inv_gaussian_family <- function(){
    g <- function(x){
        tmp <- qnorm(1 - x)
        pmax(pmin(tmp, qnorm(1 - 10^-15)), qnorm(10^-15))
    }
    ginv <- function(x){
        tmp <- 1 - pnorm(x)
        pmax(pmin(tmp, 1 - 10^-15), 10^-15)
    }
    eta <- function(mu){mu}
    mustar <- 0
    A <- function(mu){mu^2/2}
    name <- "inv_gaussian"
    family <- gaussian()

    exp_family(g, ginv, eta, mustar, A, name, family)
}
