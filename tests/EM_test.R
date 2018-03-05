load("../data/estrogen.RData")
pvals <- as.numeric(estrogen$pvals)
x <- data.frame(x = as.numeric(estrogen$ord))
n <- length(pvals)
s <- rep(0.45, n)
dist <- beta_family()
piargs <- muargs <- list(formula = ~ ns(x, df = 8))

source("source_all.R")

init_obj <- init_mix(x, pvals, s, dist, "glm",
                     piargs = piargs, muargs = muargs)
summary(init_obj$pix)
summary(init_obj$mux)
Estep_obj <- Estep_mix(pvals, s, dist, init_obj$pix, init_obj$mux)
summary(Estep_obj$Hhat)
summary(Estep_obj$phat)
Mstep_obj <- Mstep_mix(x, Estep_obj$Hhat, Estep_obj$phat, dist,
                       "glm", piargs = piargs, muargs = muargs)
summary(Mstep_obj$pix)
summary(Mstep_obj$mux)

mod <- AdaPT(x, pvals, dist, "glm", piargs = piargs, muargs = muargs)
