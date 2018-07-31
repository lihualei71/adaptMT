#### Interpolate local FDR estimates for 2D plots
interpolate <- function(x, y, z){
    x_grid <- unique(round(as.numeric(quantile(x, seq(0.01, 0.99, 0.01), na.rm = TRUE)), 4))
    x_scale <- max(x_grid) - min(x_grid)
    y_grid <- unique(round(as.numeric(quantile(y, seq(0.01, 0.99, 0.01), na.rm = TRUE)), 4))
    y_scale <- max(y_grid) - min(y_grid)
    xy_grid <- expand.grid(x_grid, y_grid)
    z <- mapply(function(xx, yy){
        dist <- sqrt((x - xx)^2 / x_scale^2 +
                         (y - yy)^2 / y_scale^2)
        ind <- which.min(dist)
        z[ind]
    }, xy_grid[, 1], xy_grid[, 2])
    z <- matrix(z, nrow = length(x_grid))
    return(list(x = x_grid, y = y_grid, z = z))
}

#' Plotting Functions for AdaPT with 2D Covariates
#'
#' Plotting the outputs of \code{adapt} when \code{x} is 2-dimensional, including threshold curves and level curves of local FDR.
#'
#' The breaks in the legend of \code{plot_2d_thresh} correspond to the maximum, the 95% quantile, the 90% quantile, the 85% quantile, the 75% quantile and the minimum of the thresholds.
#'
#' \code{plot_2d_lfdr} gives the contour plot of local FDR estimates when all p-values are equal to \code{targetp}. It is recommended to run \code{plot_2d_lfdr} for multiple \code{targetp}'s ranging from {0.001, 0.005, 0.01, 0.05}.
#'
#' @param obj an 'adapt' object
#' @param x covariates (i.e. side-information). Should be compatible to \code{models} and 2-dimensional.
#' @param pvals a vector of values in [0, 1]. P-values
#' @param alpha a positive scalar in (0, 1). Target FDR level
#' @param title a string. Title of the figure
#' @param targetp a real in (0, 1). See Details
#' @param xlab,ylab a string. Label of x/y-axis
#' @param keyaxes a list of arguments passed into axis. The graphical setting for the legend bar. An empty list by default
#' @param ... other arguments passed to \code{\link[graphics]{par}}
#'
#' @name plot_2d
#'
#' @examples
#' \donttest{
#' # Generate a 2-dim x
#' n <- 400
#' x1 <- x2 <- seq(-100, 100, length.out = 20)
#' x <- expand.grid(x1, x2)
#' colnames(x) <- c("x1", "x2")
#'
#' # Generate p-values (one-sided z test)
#' # Set all hypotheses in the central circle with radius 30 to be
#' # non-nulls. For non-nulls, z~N(2,1) and for nulls, z~N(0,1).
#' H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
#' mu <- ifelse(H0, 2, 0)
#' set.seed(0)
#' zvals <- rnorm(n) + mu
#' pvals <- 1 - pnorm(zvals)
#'
#' # Run adapt_gam with a 2d spline basis
#' library("mgcv")
#' formula <- "s(x1, x2)"
#' dist <- beta_family()
#' res <- adapt_gam(x = x, pvals = pvals, pi_formulas = formula,
#'                  mu_formulas = formula, dist = dist, nfits = 5)
#'
#' # Plots
#' plot_2d_thresh(res, x, pvals, 0.3, "P-value Thresholds (alpha = 0.3)")
#' plot_2d_lfdr(res, x, pvals, 0.3, "Local FDR Estimates (alpha = 0.3, p = 0.01)", 0.01)
#' }
NULL

#' @rdname plot_2d
#'
#' @export
plot_2d_thresh <- function(obj, x, pvals,
                           alpha, title,
                           xlab = NULL, ylab = NULL,
                           keyaxes = list(),
                           ...){
    if (!"adapt" %in% class(obj)){
        stop("obj is not an 'adapt' object.")
    }

    xnames <- colnames(x)
    xlab <- ifelse(is.null(xlab), xnames[1], xlab)
    ylab <- ifelse(is.null(ylab), xnames[2], ylab)
    if (ncol(x) != 2){
        stop("x must be 2-dimensional.")
    }
    x1 <- as.numeric(x[,1])
    x2 <- as.numeric(x[,2])
    n <- length(pvals)
    dist <- obj$dist
    alphas <- obj$alphas
    if (alpha < min(alphas)){
        s <- rep(0, n)
        nrejs <- 0
    } else {
        ind <- max(which(alphas <= alpha))
        s <- obj$s[, ind]
        nrejs <- obj$nrejs[ind]
    }

    par(...)

    if (max(s) == 0){
        warning("All thresholds are zero. Return no plot")
        return()
    }
    s_levels <- unique(round(c(0, quantile(s, c(0.75, 0.8, 0.85, 0.9, 0.95)), max(s)), 3))
    colors <- c("white", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#084594")[1:length(s_levels)]
    plot_data <- interpolate(x1, x2, s)
    keyaxes_params <- c(list(side = 4, at = s_levels),
                        keyaxes)
    filled.contour(plot_data$x, plot_data$y, plot_data$z,
                   levels = s_levels,
                   col = colors, xlab = xlab, ylab = ylab,
                   main = title,
                   key.axes = do.call(axis, keyaxes_params))
}

#' @rdname plot_2d
#'
#' @export
plot_2d_lfdr <- function(obj, x, pvals,
                         alpha, title, targetp,
                         xlab = NULL, ylab = NULL,
                         keyaxes = list(),
                         ...){
    if (!"adapt" %in% class(obj)){
        stop("\'obj\' is not an \'adapt\' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }

    xnames <- colnames(x)
    xlab <- ifelse(is.null(xlab), xnames[1], xlab)
    ylab <- ifelse(is.null(ylab), xnames[2], ylab)
    if (ncol(x) != 2){
      stop("x must be 2-dimensional.")
    }
    x1 <- as.numeric(x[,1])
    x2 <- as.numeric(x[,2])
    n <- length(pvals)
    dist <- obj$dist
    params_alphas <- sapply(obj$params, function(x){x$alpha})
    if (alpha > max(params_alphas)){
        ind <- 1
    } else {
        ind <- min(which(params_alphas >= alpha))
    }
    params <- obj$params[[ind]]
    lfdr <- compute_lfdr_mix(targetp, dist = dist, params = params)

    par(...)

    lfdr_levels <- c(0, 0.05, 0.1, 0.2, 0.5, 1)
    colors <- c("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9", "white")
    plot_data <- interpolate(x1, x2, lfdr)
    keyaxes_params <- c(list(side = 4, at = lfdr_levels),
                        keyaxes)
    filled.contour(plot_data$x, plot_data$y, plot_data$z,
                   levels = lfdr_levels,
                   col = colors, xlab = xlab, ylab = ylab,
                   main = title,
                   key.axes = do.call(axis, keyaxes_params))
}
