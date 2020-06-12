#'  Fitted probabilities from \code{brm} fits
#' 
#' @description  Calculate fitted probabilities from a fitted binary regression model object.
#' 
#' @param object    A fitted object from function \code{brm}.
#' 
#' @param va.new    An optional covariate matrix to make predictions with. If omitted, the original matrix va is used.
#' 
#' @param vb.new    An optional covariate matrix to make predictions with. If vb.new is omitted but va.new is not, then vb.new is set to be equal to va.new. If both vb.new and va.new are omitted, then the original matrix vb is used.
#' 
#' @param x.new     An optional vector of x.
#' 
#' @param ...    affecting the predictions produced.
#' 
#' @return If x.new is omitted, a matrix consisting of fitted probabilities for p0 = P(y=1|x=0,va,vb) and p1 = P(y=1|x=1,va,vb).
#' 
#' If x.new is supplied, a vector consisting of fitted probabilities px = P(y=1|x=x.new,va,vb).
#' 
#' @export


predict.brm = function(object, x.new = NULL, va.new = NULL, vb.new = NULL, ...) {
    
    va = object$va
    vb = object$vb
    
    if(is.null(vb.new)){
        if(is.null(va.new)){
            vb.new = vb
        }else{
            vb.new = va.new
        }
    }
    if(is.null(va.new)) va.new = va
    
    n  = nrow(va.new)
    pa = ncol(va.new)
    pb = ncol(vb.new)
    alpha.est = object$point.est[1:pa]
    beta.est  = object$point.est[(pa+1):(pa+pb)]
    
    linear.predictors = cbind(va.new %*% alpha.est, vb.new %*% beta.est)
    if(object$param=="RR") 
        p0p1 = getProbRR(linear.predictors)
    if(object$param=="RD") 
        p0p1 = getProbRD(linear.predictors)
    if(object$param=="OR"){
        p0 = expit(linear.predictors[,2])
        or = exp(linear.predictors[,1])
        odds1 = or * (p0 / (1-p0))
        p1 = odds1 / (1+odds1)
        p0p1 = cbind(p0, p1)
    }
    colnames(p0p1) = c("p0", "p1")
    
    if(!is.null(x.new)){
        px = rep(NA, n)     
        px[x.new == 0] = p0p1[x.new == 0, 1]
        px[x.new == 1] = p0p1[x.new == 1, 2]
        return(px)
    }else{
        return(p0p1)
    }
    
} 






