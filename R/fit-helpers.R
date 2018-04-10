#---------------------------------------------------------------
# Helpers for fitting, including M-steps and initialization
#---------------------------------------------------------------
fit_pi <- function(fun, args, type = c("Mstep", "init")){
    type <- type[1]
    if (type == "Mstep"){
        fun_name <- "\'pifun\'"
    } else if (type == "init"){
        fun_name <- "\'pifun_init\'"
    }
    
    if (is.null(args)) {
        stop(paste0(fun_name, " has irregular input types. Replace it with another function or writing a wrapper of ", fun_name, " with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)"))
    }

    fit <- do.call(fun, args)
    if (is.null(fit$fitv)){
        stop(paste0(fun_name, " does not output \'fitv\'. Replace it with another function or change the name for fitted value to \'fitv\'"))
    }

    pix <- as.numeric(fit$fitv)
    if (any(is.na(pix))){
        stop("Initialization of \'pix\' has NAs")
    }
    pix <- pminmax(pix, 0, 1)

    return(list(fitv = pix, info = fit$info))
}

fit_mu <- function(fun, args, dist, type = c("Mstep", "init")){
    type <- type[1]
    if (type == "Mstep"){
        fun_name <- "\'pifun\'"
    } else if (type == "init"){
        fun_name <- "\'pifun_init\'"
    }

    if (is.null(args)) {
        stop(paste0(fun_name, " has irregular input types. Replace it with another function or writing a wrapper of ", fun_name, " with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)"))
    }

    fit <- do.call(fun, args)
    if (is.null(fit$fitv)){
        stop(paste0(fun_name, " does not output \'fitv\'. Replace it with another function or change the name for fitted value to \'fitv\'"))
    }

    mux <- as.numeric(fit$fitv)
    if (any(is.na(mux))){
        stop("Initialization of \'mux\' has NAs")
    }
    if (dist$family$family == "Gamma"){
        mux <- pmax(mux, 1)
    } else if (dist$family$family == "gaussian"){
        mux <- pmax(mux, 0)
    }

    return(list(fitv = mux, info = fit$info))
}
