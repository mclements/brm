
ValidCheck = function(param, y, x, va, vb, vc, weights, subset, est.method, 
                      optimal, max.step, thres, alpha.start, beta.start) {
    if (!is.character(param)) 
        stop("Parameter must be a character")
    if (!(param %in% c("RD", "RR", "OR"))) 
        stop("Parameter can only take RR, RD or OR")
    
    if (sum(is.na(y)) + sum(is.na(x)) + sum(is.na(va)) + sum(is.na(vb)) + 
        sum(is.na(vc)) + sum(is.na(weights)) > 0) 
        warning("Observations with missing values will be removed.")
    if (!(all(y %in% c(0, 1)))) 
        stop("y values must be either 0 or 1.")
    if (!(all(x %in% c(0, 1)))) 
        stop("x values must be either 0 or 1.")
    if (!identical(length(y), length(x), dim(va)[1], dim(vb)[1], dim(vc)[1])) 
        stop("y, x and v must have the same length (dimension)")
    
    if(!is.numeric(weights))
        stop("weights must either be NULL or take numerical values")
    if(!is.numeric(subset))
        stop("subset must either be NULL or take numerical values")
    if (!(est.method %in% c("MLE", "DR"))) 
        stop("Must use MLE or DR for estimation")
    if (!is.logical(optimal)) 
        stop("optimal must be a logical variable")
    if(!is.numeric(max.step) & !is.null(max.step))
        stop("max.step must be a number")
    if(!is.numeric(thres))
        stop("thres must be a number")
    if(!is.null(alpha.start) & length(alpha.start) != dim(va)[2])
        stop("length of alpha.start must match the dimension of va")
    if(!is.null(beta.start) & length(beta.start) != dim(vb)[2])
        stop("length of beta.start must match the dimension of vb")
    
} 
