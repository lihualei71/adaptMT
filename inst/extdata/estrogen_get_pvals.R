#----------------------------------------------------------------
# This script for obtaining 'estrogen' dataset in 'adaptMT' package.
# The code is from Ang Li and Rina Foygel Barber (2016) 'Multiple testing with the structure adaptive Benjamini-Hochberg algorithm' We take the first part of the source code from
# 'http://www.stat.uchicago.edu/~rina/sabha/gene_drug_data_example.R' and
# modified it for our purposes. 
#----------------------------------------------------------------

#----------------------------------------------------------------
# download the gene/drug data
#----------------------------------------------------------------

if (!require(GEOquery)) {
    install.packages("GEOquery")
}
library("GEOquery")
gds2324 <- getGEO("GDS2324")
eset <- GDS2eSet(gds2324, do.log2 = TRUE)
Data <- exprs(eset)

#----------------------------------------------------------------
# calculate the p-values with different orderings
#----------------------------------------------------------------

#' code for processing this data set & for accumulation tests
#' the code in this function is from the paper:
#' Ang Li and Rina Foygel Barber, 'Accumulation tests for FDR control in ordered hypothesis testing' (2015). Available from http://arxiv.org/abs/1505.07352

#' the following code is copied (with some edits) from:
#' http://www.stat.uchicago.edu/~rina/accumulationtests/gene_dosage_experiment.R

#' Two-sample t-tests on a data matrix
#' This function inputs a data matrix X with n columns. The columns indexed by g_1 and g_2 belong to sample 1 and sample 2, respectively (e.g. cases and controls). For each row of X, the function runs a two-sample t-test for that row.
ttest_mat <- function(X, g1, g2) {
    ## g1 & g2 give columns for groups 1 & 2
    n1 <- length(g1)
    n2 <- length(g2)
    means1 <- rowMeans(X[, g1])
    means2 <- rowMeans(X[, g2])
    vars1 <- rowSums((X[, g1] - means1 %*% t(rep(1, n1)))^2)/(n1 - 1)
    vars2 <- rowSums((X[, g2] - means2 %*% t(rep(1, n2)))^2)/(n2 - 1)
    sds_diff <- sqrt(vars1 / n1 + vars2 / n2)
    tstats <- (means1 - means2) / sds_diff
    dfs <- (vars1 / n1 + vars2 / n2)^2 / ((vars1 / n1)^2 / (n1 - 1) + (vars2 / n2)^2 / (n2 - 1))
    pvals <- 2 * (1 - pt(abs(tstats), dfs))
    output <- list()
    output$tstats <- tstats
    output$dfs <- dfs
    output$pvals <- pvals
    output
}

#' The next function is the same, except that instead of performing two-sided t-tests, we perform a one-sided t-test for each row of X. The signs s_i specify, for the i-th row of X, whether we are testing for a positive effect or a negative effect.
signed_ttest_mat <- function(X, g1, g2, s) {
    n1 <- length(g1)
    n2 <- length(g2)
    means1 <- rowMeans(X[, g1])
    means2 <- rowMeans(X[, g2])
    vars1 <- rowSums((X[, g1] - means1 %*% t(rep(1, n1)))^2)/(n1 - 1)
    vars2 <- rowSums((X[, g2] - means2 %*% t(rep(1, n2)))^2)/(n2 - 1)
    sds_diff <- sqrt(vars1 / n1 + vars2 / n2)
    tstats <- s * (means1 - means2)/sds_diff
    dfs <- (vars1 / n1 + vars2 / n2)^2 / ((vars1 / n1)^2 / (n1 - 1) + (vars2 / n2)^2 / (n2 - 1))
    pvals <- (1 - pt(tstats, dfs))
    output <- list()
    output$tstats <- tstats
    output$dfs <- dfs
    output$pvals <- pvals
    output
}

#' Each row of 'Data' is a gene (total: n=22283 genes). 
#' The 25 columns of 'Data' correspond to 5 trials each at 5 different dosages: columns 1-5 are zero dose (control group), followed by columns 6-10 at the lowest dose, etc. The entries of 'Data' are log-transformed gene expression levels.
gene_drug_get_pvals <- function(Data, highdose = 21:25, lowdose = 6:10, control = 1:5) {
    n <- dim(Data)[1]
    
    ## Computing p-values Next, we will use the highest-dosage data to produce an
    ## ordering on the low-dosage vs control-group p-values, resulting in an ordered
    ## sequence of p-values which we will use for the accumulation tests.  First, for
    ## each gene, produce a pval for highest dose compared to mean of lowest dose +
    ## control. Record also the sign of the effect (increased or reduced gene
    ## expression at the highest dose, compared to pooled low-dose and control-group
    ## data).
    ttest_highdose <- ttest_mat(Data, highdose, c(lowdose, control))
    pvals_highdose <- ttest_highdose$pvals
    test_signs <- sign(ttest_highdose$tstats)
    
    pvals_highdose_small <- ttest_mat(Data, highdose[1:2], c(lowdose, control))$pvals
    set.seed(1)
    
    ## Next, for each gene we will perform a one-sided t-test (using the sign of the
    ## estimated high-dose effect), and then use a permutation test to get a final
    ## p-value for this gene. These p-values will then be reordered according to the
    ## high-dose results above.
    signed_pvals_lowdose_ttest <- signed_ttest_mat(Data, lowdose, control, test_signs)$pvals
    
    nhigh <- length(highdose)
    nlow <- length(lowdose)
    ncontrol <- length(control)
    nsmall <- nlow + ncontrol
    smalldose <- c(lowdose, control)
    signed_pvals_lowdose_ttest_permuted <- matrix(0, choose(nsmall, nlow), n)
    nchoosek <- combn(smalldose, nlow)
    for (i in 1:choose(nsmall, nlow)) {
        permord <- c(nchoosek[, i], setdiff(smalldose, nchoosek[, i]))
        temp_data <- Data
        temp_data[, smalldose] <- temp_data[, permord]
        signed_pvals_lowdose_ttest_permuted[i, ] <- signed_ttest_mat(temp_data, lowdose, 
            control, test_signs)$pvals
    }
    signed_pvals_permutation_test <- colSums(abs(signed_pvals_lowdose_ttest_permuted) - 
        rep(1, choose(nsmall, nlow)) %*% t(abs(signed_pvals_lowdose_ttest)) <= 0)/choose(nsmall, 
        nlow)
    n <- nrow(Data)
    pvals <- signed_pvals_permutation_test - 1/choose(nsmall, nlow) * runif(n)
    ## multiplying to avoid pvalues exactly equal to 1 due to discretization
    ord <- rank(abs(pvals_highdose))
    ord_small <- rank(abs(pvals_highdose_small))
    
    output <- data.frame(pvals = pvals, ord = ord, ord_small = ord_small)
    return(output)
}

estrogen <- gene_drug_get_pvals(Data)
save(file = "../../data/estrogen.RData", estrogen)
