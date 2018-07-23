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
                 verbose = list(print = FALSE, fit = FALSE, ms = FALSE),
                 s0 = rep(0.5, length(pvals)))

## Begin Tests
test_that("length of 'order' should match the number of hypotheses (#1)", {
    n <- length(pvals)
    expect_equal(length(res$order), n)
    expect_equal(length(unique(res$order)), n)
})

test_that("'order' should be consistent with thresholds 's'", {
    for (i in 1:100){
        mask <- which(pvals <= res$s[, i] | pvals >= 1 - res$s[, i])
        mask2 <- sort(tail(res$order, length(mask)))
        diffs <- setdiff(mask, mask2)
        expect_equal(length(diffs), 0)
    }
})

test_that("'nrejs' should be consistent with thresholds 's'", {
    nrejs2 <- apply(res$s, 2, function(x){
        sum(pvals <= x)
    })
    expect_equal(res$nrejs, nrejs2)
})

test_that("'nrejs' should be consistent with 'rejs'", {
    nrejs3 <- sapply(res$rejs, length)
    expect_equal(res$nrejs, nrejs3)
})

test_that("'s' is decreasing", {
    ncol <- ncol(res$s)
    s_diff <- res$s[, 2:ncol] - res$s[, 1:(ncol - 1)] + 1e-10
    expect_equal(all(s_diff >= 0), TRUE)
})

test_that("'qvals' is consistent with 'rejs", {
    for (i in 1:100){
        rejs1 <- res$rejs[[i]]
        rejs2 <- which(res$qvals <= res$alphas[i])
        expect_equal(rejs1, rejs2)
    }
})
