#---------------------------------------------------------------
# Helpers
#---------------------------------------------------------------

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

complete_pkg <- function(formula){
    formula <- as.character(formula)
    formula <- tail(formula, 1)
    formula <- tail(strsplit(formula, "~")[[1]], 1)    
    formula <- paste0(" ", formula)    
    if (grepl("ns\\(", formula)){
        if (!requireNamespace("splines", quietly = TRUE)){
            stop("package \'splines\' not found. Please install.")
        }
        formula <- gsub("([^:])ns\\(", "\\1splines::ns\\(", formula)
    }
    if (grepl("[^a-z]s\\(", formula)){
        if (!requireNamespace("mgcv", quietly = TRUE)){
            stop("package \'mgcv\' not found. Please install.")
        }
        formula <- gsub("([^:a-z])s\\(", "\\1mgcv::s\\(", formula)
    }
    return(formula)
}
    

complete_formula <- function(formula, response_name){
    if (is.null(formula)){
        stop("No formula is found. Please specify a formula ")
    }
    formula <- as.character(formula)
    formula <- tail(formula, 1)
    formula <- tail(strsplit(formula, "~")[[1]], 1)
    formula <- paste0(" ", formula)
    ## completed_formula <- as.formula(
    ##     paste(response_name, "~", formula),
    ##     env = environment(args$formula))
    completed_formula <- paste0(response_name, " ~", formula)

    return(completed_formula)
}

complete_args <- function(x, response, fun,
                          args = NULL,
                          weights = NULL,
                          force_integer = FALSE){
    input_type <- func_input_type(fun)
    if (!input_type %in% c("formula", "xy", "Xy")){
        stop("Wrong input type.")
    }
    
    response_name <- find_newname(colnames(x))

    if (input_type == "formula"){
        if (is.null(args) || !"formula" %in% names(args)){
            stop("Formula is not found. Please specify a formula for the fitting function.")
        }
        data <- cbind(data.frame(response), x)        
        colnames(data)[1] <- response_name
        args$formula <-  complete_formula(args$formula, response_name)
        data_args <- c(list(data = data), args)
    } else if (input_type == "xy"){
        data_args <- c(
            list(x = x, y = response),
            args)
    } else if (input_type == "Xy"){
        data_args <- c(
            list(X = x, y = response),
            args)
    } 

    data_args <- c(data_args, list(weights = weights))

    return(data_args)
}

complete_model <- function(model, dist){
    if (is.null(model$algo)){
        switch(model$name,
               "glm" = gen_adapt_model_glm(
                   dist, model$args$piargs, model$args$muargs
                   ),
               "gam" = gen_adapt_model_gam(
                   dist, model$args$piargs, model$args$muargs
                   ),
               "glmnet" = gen_adapt_model_glmnet(
                   dist, model$args$piargs, model$args$muargs
                   ),
               stop("\'model$name\' not found in the library")
               )
    } else {
        model
    }
}
