logit <- function(x){
    log(x / (1 - x))
}

inv.logit <- function(x){
    exp(x) / (1 + exp(x))
}

#' Function to calculate FDPhat
#' 
#' @param pvals: p-values
#' @param s: threshold curve
#' @return 
#'    fdphat: estimated FDP
fdp.hat <- function(pvals, s) {
    A <- sum(pvals > 1 - s)
    R <- sum(pvals <= s)
    return((1 + A) / max(R, 1))
}

func_input_type <- function(fun){
    argnames <- formalArgs(fun)
    if ("formula" %in% argnames){
        return("formula")
    } else if ("x" %in% argnames){
        return("xy")
    } else if ("X" %in% argnames){
        return("Xy")
    }
}

find_newname <- function(names.vec){
    name <- "aaa"
    while (name %in% names.vec){
        name <- paste0(name, "a")
    }
    return(name)
}

complete_formulas <- function(params, response.name){
    if (!is.null(params[["formula"]])){
        completed.formula <- as.formula(
            paste(response.name, "~", params[["formula"]]))
        params[["formula"]] <- completed.formula
    }
    if (!is.null(params[["alter.formulas"]])){
        completed.alter.formulas <-
            sapply(params[["alter.formulas"]], function(formula){
                as.formula(paste(response.name, "~", formula))
            })
        params[["alter.formulas"]] <- completed.alter.formulas
    }
    return(params)
}
