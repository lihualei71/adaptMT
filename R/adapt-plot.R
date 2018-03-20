adapt_1d_plot <- function(obj, alpha, title,
                          data = NULL,
                          xlab = "x",
                          disp.pmax = 0.2,
                          num.xbreaks = 3,
                          rand.seed.purturb = NA,
                          legend.position = "topright"){
    if (!"adapt" %in% class(obj)){
        stop("\"obj\" is not of class \"adapt\".")
    }

    if (is.null(data)){
        data <- obj$data
    }
    x <- as.numeric(data[["x"]][,1])
    pvals <- data[["pvals"]]
    dist <- obj$info$dist
    alphas <- obj$alphas
    ind <- which.min(abs(alphas - alpha))
    s <- obj$s[, ind]
    pix <- obj$params[[ind]]$pix
    mux <- obj$params[[ind]]$mux
    n <- length(s)
    
    par(mfrow = c(2, 1))
    par(mar=c(4.1, 4.1, 2, 0.15),
        cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)

    ## Top panel on threshold curves
    what.type <- ifelse(pvals < s, 1, ifelse(pvals > 1 - s, 2, 3))
    plot(x, (1:n * disp.pmax) / n, type = "n", pch = ".",
         xaxs = "i", yaxs = "i", ylab = "p-value", xlab = "",
         col = c("red", "blue", "black")[what.type], yaxt = "n",
         main = paste("Rejection threshold (", title, ")",
             sep = ""))
    axis(2, at = seq(0, disp.pmax, length.out = num.xbreaks),
         labels = seq(0, disp.pmax, length.out = num.xbreaks))
    hhi <- hlo <- .5
    vlo <- -.01; vhi <- 1.01
    polygon(x = c(x, rev(x)), y = c(s, rep(vlo, n)),
            col = "#FFDDDD", border = "red")
    polygon(x = c(x, rev(x)), y = c(1 - s, rep(vhi, n)),
            col = "light blue", border = "blue")
    if(!is.na(rand.seed.purturb)) {
        set.seed(rand.seed.purturb)
        points(x, pvals + 0.001 * rnorm(n),
               pch = ".",
               col = c("red", "blue", "black")[what.type])
    }
    points(x, pvals, pch = ".",
           col = c("red", "blue", "black")[what.type])
    box()

    ## Bottom panel on local fdr level curves.
    ind <- order(x)
    x <- x[ind]
    mux <- mux[ind]
    pix <- pix[ind]    

    x.grid <- c(seq(1, n, by = 100))
    if (x.grid[length(x.grid)] < n){
        x.grid <- c(x.grid, n)
    }
    p.grid <- exp(seq(-15, log(disp.pmax), length.out=100))
    xp.grid <- expand.grid(x = x.grid, p = p.grid)
    locfdr.vals <- lfdr.mix(p = xp.grid$p, dist = dist,
                            params = list(mux = mux[xp.grid$x],
                                          pix = pix[xp.grid$x]))
    locfdr.mat <- matrix(locfdr.vals, nrow = length(x.grid))
    plot(0, 0, xlim = c(min(x), max(x)), ylim = c(0, disp.pmax),
         type = "n", xaxs = "i", yaxs = "i", yaxt = "n", 
         main = paste("Estimated local FDR (", title, ")", sep=""),
         ylab = "p-value", xlab = xlab)
    axis(2, at = seq(0, disp.pmax, length.out = num.xbreaks),
         labels = seq(0, disp.pmax, length.out = num.xbreaks))
    colors <- c(rev(brewer.pal(n = 4, name = "Reds")), "white")
    .filled.contour(x = x[x.grid] + 1e-10 * (1:length(x.grid)),
                    y = p.grid, z = locfdr.mat,
                    levels = c(0, 0.1, 0.2, 0.3, 0.5, 1),
                    col = colors)
    legend(legend.position, col = "black", fill = rev(colors),
           legend = rev(c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.5","0.5-1")),
           cex = 1.3, bty = "n")
    
}
