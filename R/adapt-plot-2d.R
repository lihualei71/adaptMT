#### Interpolate local FDR estimates for 2-d plots
interpolate <- function(x, y, z){
    x_grid <- unique(round(as.numeric(quantile(x, seq(0.01, 0.99, 0.01), na.rm = TRUE)), 4))
    x_scale <- max(x_grid) - min(x_grid)
    y_grid <- unique(as.numeric(quantile(y, seq(0.01, 0.99, 0.01), na.rm = TRUE)))
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

#' Plotting Functions for AdaPT with 2d Covariates
#'
#' Plotting the outputs of \code{adapt} when \code{x} is 2-dimensional, including threshold curves and level curves of local FDR.
#'
#' The breaks in the legend of \code{plot_2d_thresh} correspond to the maximum, the 95% quantile, the 90% quantile, the 85% quantile, the 75% quantile and the minimum of the thresholds.
#'
#' \code{plot_2d_lfdr} gives the contour plot of local FDR estimates when all p-values are equal to \code{targetp}. It is recommended to run \code{plot_2d_lfdr} for multiple \code{targetp}'s ranging from {0.001, 0.005, 0.01, 0.05}. 
#' 
#' @param obj an 'adapt' object
#' @param alpha a positive scalar in (0, 1). Target FDR level
#' @param title a string. Title of the figure
#' @param targetp a real in (0, 1). See Details
#' @param data a list in the form of list(x = , pvals = ). NULL if obj$data is not NULL
#' @param xlab,ylab a string. Label of x/y-axis 
#' @param ... other arguments passed to \code{\link[graphics]{par}}
#'
#' @name plot_2d
#'
NULL

#' @rdname plot_2d
#'
#' @export
plot_2d_thresh <- function(obj, alpha, title,
                           data = NULL,
                           xlab = NULL, ylab = NULL,
                           keyaxes = list(),
                           ...){
    if (!"adapt" %in% class(obj)){
        stop("obj is not an 'adapt' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }
    xnames <- colnames(data[["x"]])
    xlab <- ifelse(is.null(xlab), xnames[1], xlab)
    ylab <- ifelse(is.null(ylab), xnames[2], ylab)
    x1 <- as.numeric(data[["x"]][,1])
    x2 <- as.numeric(data[["x"]][,2])
    pvals <- data[["pvals"]]
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

    s_levels <- round(c(0, quantile(s, c(0.75, 0.8, 0.85, 0.9, 0.95)), max(s)), 3)
    colors <- c("white", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#084594")
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
plot_2d_lfdr <- function(obj, alpha, title, targetp,
                         data = NULL,
                         xlab = NULL, ylab = NULL,
                         keyaxes = list(),
                         ...){
    if (!"adapt" %in% class(obj)){
        stop("\'obj\' is not an \'adapt\' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }
    xnames <- colnames(data[["x"]])
    xlab <- ifelse(is.null(xlab), xnames[1], xlab)
    ylab <- ifelse(is.null(ylab), xnames[2], ylab)
    x1 <- as.numeric(data[["x"]][,1])
    x2 <- as.numeric(data[["x"]][,2])
        pvals <- data[["pvals"]]
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
