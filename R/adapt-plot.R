#' Plotting
#' 
#' @export
plot_thresh_1d <- function(obj, alpha, title,
                           data = NULL,
                           xlab = "x",
                           disp_ymax = 0.2,
                           num_xbreaks = 3,
                           rand_seed_perturb = NA,
                           ...){
    if (!"adapt" %in% class(obj)){
        stop("obj is not an 'adapt' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }
    x <- as.numeric(data[["x"]][,1])
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

    x_ord <- order(x)
    pvals <- pvals[x_ord]
    s <- s[x_ord]
    x <- x[x_ord]

    par(...)
    ## par(mfrow = c(2, 1), ...)

    ## Top panel on threshold curves
    what_type <- ifelse(pvals < s, 1, ifelse(pvals > 1 - s, 2, 3))
    plot(x, (1:n * disp_ymax) / n, type = "n", pch = ".",
         xaxs = "i", yaxs = "i", ylab = "p-values", xlab = xlab,
         col = c("red", "blue", "black")[what_type], yaxt = "n",
         main = title)
    axis(2, at = seq(0, disp_ymax, length.out = num_xbreaks),
         labels = seq(0, disp_ymax, length.out = num_xbreaks))
    hhi <- hlo <- .5
    vlo <- -.01
    vhi <- 1.01
    polygon(x = c(x, rev(x)), y = c(s, rep(vlo, n)),
            col = "#FFDDDD", border = "red")
    polygon(x = c(x, rev(x)), y = c(1 - s, rep(vhi, n)),
            col = "light blue", border = "blue")
    if(!is.na(rand_seed_perturb)) {
        set.seed(rand_seed_perturb)
        points(x, pvals + 0.001 * rnorm(n),
               pch = ".",
               col = c("red", "blue", "black")[what_type])
    }
    points(x, pvals, pch = ".",
           col = c("red", "blue", "black")[what_type])
    box()
}

#' @export
plot_lfdr_1d <- function(obj, alpha, title,
                         data = NULL,
                         xlab = "x",
                         disp_ymax = 0.2,
                         num_xbreaks = 3,
                         legend_pos = "topright",
                         ...){
    if (!"adapt" %in% class(obj)){
        stop("obj is not an 'adapt' object.")
    }

    if (is.null(data)){
        data <- obj$data
    }
    x <- as.numeric(data[["x"]][,1])
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

    x_ord <- order(x)
    pvals <- pvals[x_ord]
    pix <- pix[x_ord]
    mux <- mux[x_ord]
    x <- x[x_ord]

    par(...)
    
    ## Bottom panel on local fdr level curves.
    x_grid <- c(seq(1, n, by = 100))
    if (x_grid[length(x_grid)] < n){
        x_grid <- c(x_grid, n)
    }
    p_grid <- exp(seq(-15, log(disp_ymax), length.out=100))
    xp_grid <- expand.grid(x = x_grid, p = p_grid)
    locfdr_vals <-
        compute_lfdr_mix(pvals = xp_grid$p, dist = dist,
                         params = list(mux = mux[xp_grid$x],
                             pix = pix[xp_grid$x]))
    locfdr_mat <- matrix(locfdr_vals, nrow = length(x_grid))
    plot(0, 0, xlim = c(min(x), max(x)), ylim = c(0, disp_ymax),
         type = "n", xaxs = "i", yaxs = "i", yaxt = "n", 
         main = title, ylab = "p-value", xlab = xlab)
    axis(2, at = seq(0, disp_ymax, length.out = num_xbreaks),
         labels = seq(0, disp_ymax, length.out = num_xbreaks))
    colors <- c("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9", "white")
    .filled.contour(x = x[x_grid] + 1e-10 * (1:length(x_grid)),
                    y = p_grid, z = locfdr_mat,
                    levels = c(0, 0.1, 0.2, 0.3, 0.5, 1),
                    col = colors)
    legend(legend_pos, col = "black", fill = rev(colors),
           legend = rev(c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.5","0.5-1")),
           bty = "n")
    
}
