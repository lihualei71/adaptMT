#' @export

adapt_glm <- function(x, pvals, pi_formulas, mu_formulas,
                      dist = beta_family(),
                      s0 = rep(0.45, length(pvals)),
                      alphas = seq(0.01, 1, 0.01),
                      ...){

    models <- lapply(1:length(pi_formulas), function(i){
        piargs <- list(formula = pi_formulas[[i]])
        muargs <- list(formula = mu_formulas[[i]])
        gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
    })

    adapt(x, pvals, models, dist, s0, alphas, ...)
}
