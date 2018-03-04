pminmax <- function(x, low, up){
    pmin(pmax(x, low), up)
}

logit <- function(x){
    log(x / (1 - x))
}

inv_logit <- function(x){
    exp(x) / (1 + exp(x))
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

find_newname <- function(names_vec){
    name <- "aaa"
    while (name %in% names_vec){
        name <- paste0(name, "a")
    }
    return(name)
}

complete_formulas <- function(args, response_name){
    if (!is.null(args[["formula"]])){
        formula <- as.character(args[["formula"]])
        formula <- tail(formula, 1)
        formula <- tail(strsplit(formula, "~")[[1]], 1)
        completed_formula <- as.formula(
            paste(response_name, "~", formula))
        args[["formula"]] <- completed_formula
    } else {
        stop("No formula is found. Please specify a formula ")
    }
    if (!is.null(args[["alter.formulas"]])){
        completed_alter_formulas <-
            sapply(args[["alter.formulas"]], function(formula){
                formula <- as.character(args[["formula"]])
                formula <- tail(formula, 1)
                formula <- tail(strsplit(formula, "~")[[1]], 1)
                as.formula(paste(response_name, "~", formula))
            })
        args[["alter.formulas"]] <- completed_alter_formulas
    }
    return(args)
}

complete_args <- function(x, response, fun,
                          args = NULL,
                          weights = NULL){
    input_type <- func_input_type(fun)
    response_name <- find_newname(colnames(x))

    if (input_type == "formula"){
        if (is.null(args) || !"formula" %in% names(args)){
            stop("Formula is not found. Please specify a formula for the fitting function.")
        }
        data <- cbind(data.frame(response), x)
        colnames(data)[1] <- response_name
        args <-  complete_formulas(args, response_name)
        data_args <- c(list(data = data), args)
    } else if (input_type == "xy"){
        data_args <- c(list(x = x, y = response), args)
    } else if (input_type == "Xy"){
        data_args <- c(list(X = x, y = response), args)
    } else {
        data_args <- NULL
    }

    if (!is.null(weights) && !is.null(data_args)){
        data_args$weights <- weights
    }

    return(data_args)
}
