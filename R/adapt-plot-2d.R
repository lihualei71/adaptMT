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

#' Plotting for 2d Covariate
#'
#' Plotting of the outputs of \code{adapt} when \code{x} is 2-dimensional, including threshold curves and level curves of local FDR.
#'
#' @param obj an 'adapt' object
#' @param alpha a positive scalar in (0, 1). Target FDR level
#' @param title a string. Title of the figure
#' @param data a list in the form of list(x = , pvals = ). NULL if obj$data is not NULL
#' @param xlab a string. Label of the x-axis
#' @param xlim a vector of length 2. Limits of x-axis
#' @param disp_ymax a positive scalar in (0, 1]. Maximum value displayed in the y-axis
#' @param num_yticks a positive integer. Number of ticks in the y-axis
#' @param rand_seed_perturb random seed if jitter is added. NA if no jittering is needed
#' @param legend_pos a string. Position of the legend
#' @param ... other arguments passed to \code{\link[graphics]{par}}
#'
#' @name plot_1d
#'
NULL

#' @rdname plot_1d
#'
#' @export
plot_2d_thresh <- function(obj, alpha, title,
                           data = NULL,
                           xlab = NULL, ylab = NULL,
                           ...){
    if (!"adapt" %in% class(obj)){
        stop("obj is not an 'adapt' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }
    xnames <- colnames(data[["x"]])
    xlab <- ifelse(is.null(xlab), xnames[1], xlab)
    ylab <- ifelse(is.null(xlab), xnames[2], xlab)
    x1 <- as.numeric(data[["x"]][,1])
    x2 <- as.numeric(data[["x"]][,2])
    pvals <- data[["pvals"]]
    n <- length(pvals)
    dist <- obj$dist
    alphas <- obj$alphas
    if (alpha < min(alphas)){
        s <- rep(0, n)
    } else {
        ind <- max(which(alphas <= alpha))
        s <- obj$s[, ind]
    }

    par(...)

    
}

#' @rdname plot_2d
#'
#' @export
plot_2d_lfdr <- function(obj, alpha, title, targetp,
                         data = NULL,
                         xlab = NULL, ylab = NULL,
                         ...){
    if (!"adapt" %in% class(obj)){
        stop("\'obj\' is not an \'adapt\' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }
    xnames <- colnames(data[["x"]])
    xlab <- ifelse(is.null(xlab), xnames[1], xlab)
    ylab <- ifelse(is.null(xlab), xnames[2], xlab)
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
    pix <- obj$params[[ind]]$pix
    mux <- obj$params[[ind]]$mux
    locfdr

    if (!is.null(xlim)){
        inds <- which(x >= xlim[1] & x <= xlim[2])
        x <- x[inds]
        pvals <- pvals[inds]
        pix <- pix[inds]
        mux <- mux[inds]
        n <- length(inds)
    } else {
        xlim <- c(min(x), max(x))
    }

    par(...)

    ## Bottom panel on local fdr level curves.
    x_grid <- unique(floor(seq(1, n, length.out = 100)))
    if (x_grid[length(x_grid)] < n){
        x_grid <- c(x_grid, n)
    }
    p_grid <- exp(seq(-15, log(disp_ymax), length.out=100))
    xp_grid <- expand.grid(x = x_grid, p = p_grid)
    locfdr_vals <-
        compute_lfdr_mix(pvals = xp_grid$p,
                         dist = dist,
                         params = list(
                             mux = mux[xp_grid$x],
                             pix = pix[xp_grid$x])
                         )
    locfdr_mat <- matrix(locfdr_vals, nrow = length(x_grid))
    plot(0, 0, xlim = xlim, ylim = c(0, disp_ymax),
         type = "n", xaxs = "i", yaxs = "i", yaxt = "n",
         main = title, ylab = "p-value", xlab = xlab)
    axis(2, at = seq(0, disp_ymax, length.out = num_yticks),
         labels = seq(0, disp_ymax, length.out = num_yticks))
    colors <- c("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9", "white")
    .filled.contour(x = x[x_grid] + 1e-10 * (1:length(x_grid)),
                    y = p_grid, z = locfdr_mat,
                    levels = c(0, 0.1, 0.2, 0.3, 0.5, 1),
                    col = colors)
    legend(legend_pos, col = "black", fill = rev(colors),
           legend = rev(c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.5","0.5-1")),
           bty = "n")

}
