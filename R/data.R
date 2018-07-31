#' Gene/Drug response dataset
#'
#' P-values and ordering of genes drawn from a microarray dataset, consisting of 22283 genes on breast cancer cells in response to estrogen, from NCBI Gene Expression Omnibus (GEO) through 'GEOquery' package, with index "GDS2324".
#'
#' The original dataset "GDS2324" consists of gene expression measurements for n = 22283 genes, in response to estrogen treatments in breast cancer cells for five groups of patients, with different dosage levels and 5 trials in each. The task is to identify the genes responding to a low dosage. The p-value for gene i is obtained by a one-sided permutation test which evaluates evidence for a change in gene expression level between the control group (placebo) and the low-dose group. The p-values are then ordered according to permutation t-statistics comparing the control and low-dose data, pooled, against data from a higher dosage (with genes that appear to have a strong response at higher dosages placed earlier in the list).
#'
#' Two orderings are considered: first, a stronger (more informative) ordering based on a comparison to the highest dosage; and second, a weaker (less informative) ordering based on a comparison to a medium dosage.
#'
#' The variables are as follows:
#' \itemize{
#' \item pvals. p-values
#' \item ord_high. stronger ordering
#' \item ord_mod. weaker ordering
#' }
#'
#' The R code to produce the data can be found in '/extdata/estrogen_get_pvals.R'.
#'
"estrogen"
