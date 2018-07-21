context("Test AdaPT")

## Preparation
# Load estrogen data
data(estrogen)
pvals <- as.numeric(estrogen$pvals)
x <- data.frame(x = as.numeric(estrogen$ord_high))
dist <- beta_family()

# Subsample the data for convenience
inds <- (x$x <= 5000)
pvals <- pvals[inds]
x <- x[inds,,drop = FALSE]

# Generate models for function adapt
library("splines")
formula <- "ns(x, df = 6)"

# Run adapt
res <- adapt_glm(x = x, pvals = pvals, pi_formula = formula, mu_formula = formula, nfits = 5,
                 verbose = list(print = FALSE, fit = FALSE, ms = FALSE))

## Begin Tests
test_that("length of 'order' and 'fdp' should match the number of hypotheses (#1)", {
    expect_equal(length(res$order), 5000)
    expect_equal(length(unique(res$order)), 5000)
    expect_equal(length(res$fdp), 5000)
})

test_that("'order' should be consistent with thresholds 's'", {
    for (i in 1:100){
        mask <- which(pvals <= res$s[, i] | pvals >= 1 - res$s[, i])
        mask2 <- sort(tail(res$order, length(mask)))
        expect_equal(mask, mask2)
    }
})

test_that("'nrejs' should be consistent with thresholds 's'", {
    nrejs2 <- apply(res$s, 2, function(x){
        sum(pvals <= x)
    })
    expect_equal(res$nrejs, nrejs2)
})

