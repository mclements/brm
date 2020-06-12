
#' Ancillary function for printing
#'
#' @param x a list obtained with the function 'brm'
#' 
#' @param ... additional arguments affecting the output
#' 
#' @export

print.brm = function(x, ...) {
    hid = attr(x, "hidden")
    nhid = which(!names(x) %in% hid)
    
    if (x$param == "RR") {
        cat("Parameter of interest: (conditional) relative risk;", "\n", "nuisance parameter: odds product.", 
            "\n\n", sep = "")
        cat("Target model:   log(RR) = alpha * va", "\n")
        cat("Nuisance model: log(OP) = beta * vb", "\n\n")
    }
    if (x$param == "RD") {
        cat("Parameter of interest: (conditional) risk difference;", "\n", 
            "nuisance parameter: odds product.", "\n\n", sep = "")
        cat("Target model:   log(RD) = alpha * va", "\n")
        cat("Nuisance model: log(OP) = beta * vb", "\n\n")
    }
    if (x$param == "OR") {
        cat("Parameter of interest: (conditional) odds ratio;", "\n", "nuisance parameter: baseline risk.", 
            "\n\n", sep = "")
        cat("Target model:   log(OR) = alpha * va", "\n")
        cat("Nuisance model: log(p0) = beta * vb", "\n\n")
    }
    
    
    for (i in nhid) {
        x[[i]] = round(x[[i]], 3)
    }
    
    print(x[nhid], 3)
    
    cat("See the element '$coefficients' for more information.\n")
}
