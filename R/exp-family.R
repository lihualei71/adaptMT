#===============================================================
# exp_family class
#===============================================================

#' Generate exp_family Objects for Exponential Families
#'
#' \code{exp_family} objects contain all required information in an exponential family to perform the E-step. The exponential function is encoded by
#' \deqn{h(p; \mu) = \exp\{(\eta(\mu) - \eta(\mu^{*})) g(p) - (A(\mu) - A(\mu^{*}))\}}{h(p; \eta) = exp{(\eta(\mu) - \eta(\mu*)) g(p) - (A(\mu) - A(\mu*))}}
#' where \eqn{g(p)} is an arbitrary transformation, \eqn{\mu} is the
#' \emph{mean parameter}, \eqn{\eta} is the natural parameter,
#' and \eqn{A(\mu)} is the partition function. The extra redundant
#' parameter \eqn{\mu^{*}}{\mu*} is to guarantee that \eqn{U([0, 1])}
#' belongs to the class.
#'
#' Beta family (\code{beta_family()}): modeling p-values as Beta-distributed random variables, i.e. \eqn{g(p) = -log(p)}, \eqn{\eta(\mu) = -1 / \mu}, \eqn{\mu* = 1}, \eqn{A(\mu) = log(\mu)}, name = "beta" and family = Gamma(). Beta-family is highly recommended for general problems and used as default.
#'
#' Inverse-gaussian family (\code{inv_gaussian_family()}): modeling p-values as transformed z-scores, i.e. \eqn{g(p) = \Phi^{-1}(p) (\Phi is the c.d.f. of a standard normal random variable)}, \eqn{\eta(\mu) = \mu}, \eqn{\mu* = 0}, \eqn{A(\mu) = \mu^2 / 2}, name = "inv_gaussian" and family = gaussian().
#'
#' @param g a function. An transformation of p-values
#' @param ginv a function. The inverse function of \code{g}
#' @param eta a function. The natural parameter as a function of the mean parameter \code{mu}
#' @param mustar a scalar. The mean parameter that gives \eqn{U([0, 1])}
#' @param A a function. The partition function
#' @param name a string. A name for the family. NULL by default
#' @param family an object of class "\code{\link[stats]{family}}" from \code{stats} package. The family used for model fitting in \code{\link[stats]{glm}}, \code{\link[mgcv]{gam}}, \code{\link[glmnet]{glmnet}}, etc..
#' @return an object of class "exp_family". This includes all inputs and  \code{h}, the density function.
#'
#' @export
#'
gen_exp_family <- function(g, ginv, eta, mustar, A,
                           name = NULL, family = NULL){
    h <- function(p, mu){
        ifelse(mu == mustar,
               rep(1, length(p)),
               exp(
                   g(p) * (eta(mu) - eta(mustar)) -
                       (A(mu) - A(mustar))
                   )
               )
    }
    result <- structure(
        list(h = h,
             g = g,
             ginv = ginv,
             eta = eta,
             mustar = mustar,
             A = A,
             name = name,
             family = family),
        class = "exp_family"
        )
    return(result)
}

#' @rdname gen_exp_family
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

    gen_exp_family(g, ginv, eta, mustar, A, name, family)
}

#' @rdname gen_exp_family
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

    gen_exp_family(g, ginv, eta, mustar, A, name, family)
}
