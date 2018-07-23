# adaptMT

[![Build Status](https://travis-ci.org/lihualei71/adaptMT.svg?branch=master)](https://travis-ci.org/lihualei71/adaptMT)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/lihualei71/adaptMT?branch=master&svg=true)](https://ci.appveyor.com/project/lihualei71/adaptMT)

## Overview
This package implements Adaptive P-Value Thresholding in the paper: [AdaPT: An interactive procedure for multiple testing with side information](https://arxiv.org/abs/1609.06035). It includes both a framework that allows the user to specify any algorithm to learn local FDR and a pool of convenient functions that implement specific algorithms:

- `adapt()` provides a generic framework of AdaPT permitting any learning algorithm;
- `adapt_glm()`, `adapt_gam()` and `adapt_glmnet()` provide convenient wrappers of AdaPT using Generalized Linear Models (GLM), Generalized Additive Models (GAM) and L1-penalized GLMs;

Install the adaptMT package then read vignette("adapt_demo", package = "adaptMT").

## Installation         

```
# install.packages("devtools")
devtools::install_github("lihualei71/adaptMT")
```
If one wants to access the vignette, run the following code to build the vignette. This might update other related packages and please be patient if so.

```
devtools::install_github("lihualei71/adaptMT", build_vignettes = TRUE)
```

### An Example
We illustrate the usage of adaptMT package using the example discussed in Section 5.1 of the paper [AdaPT: An interactive procedure for multiple testing with side information](https://arxiv.org/abs/1609.06035).

```
# Load package
library("adaptMT")

# Load data
data(estrogen)
pvals <- as.numeric(estrogen$pvals)
x <- data.frame(x = as.numeric(estrogen$ord_high))

# Define the exponential family for AdaPT (Section 4)
dist <- beta_family()

# Run adapt_glm
library("splines")
formulas <- paste0("ns(x, df = ", 6:10, ")")
res <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                 mu_formulas = formulas, dist = dist, nfits = 10)

# Plot the threshold curve and the level curves of local FDR
plot_1d_thresh(res, alpha = 0.1, "P-Value Thresholds")
plot_1d_lfdr(res, alpha = 0.1, "Level Curves of Local FDR Estimates")
```
