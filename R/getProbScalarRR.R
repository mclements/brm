

#' Calculate risks from log RR and log OP
#'
#' @param logrr log of relative risk
#' 
#' @param logop log of odds product
#'
#' @details The \eqn{log OP} is defined as \eqn{log OP = log[(P(y=1|x=0)/P(y=0|x=0))*(P(y=1|x=1)/P(y=0|x=1))]}. 
#' 
#' @return a vector \eqn{(P(y=1|x=0),P(y=1|x=1))}
#' 
#' @examples getProbScalarRR(0,0)
#' 
#' set.seed(0)
#' logrr = rnorm(10,0,1)
#' logop = rnorm(10,0,1)
#' probs = mapply(getProbScalarRR, logrr, logop)
#' rownames(probs) = c("P(y=1|x=0)","P(y=1|x=1)")
#' probs
#' 
#' @export

getProbScalarRR = function(logrr, logop = NA) {
    
    if(length(logrr) == 2){
        logop = logrr[2]
        logrr = logrr[1]
    }
    
    if ((logop < (-12)) || (logop > 12) || (logrr < (-12)) || (logrr > 12)) {
        # on the boundary South edge: large -ve logrr or (large -ve logop and -ve
        # logrr)
        if ((logrr < (-12)) || ((logop < (-12)) && (logrr < 0))) {
            p0 = getPrbAux(logop - logrr)
            p1 = 0
        } else {
            if ((logrr > 12) || ((logop < (-12)) && (logrr > 0))) {
                # West edge: large +ve logrr or (large -ve logop and +ve logrr)
                p0 = 0
                p1 = getPrbAux(logrr + logop)
            } else {
                # North or East edges (=logop is large +ve)
                p0 = min(exp(-logrr), 1)
                p1 = min(exp(logrr), 1)
            }
        }
    } else {
        # not on the boundary logop = 0; solving linear equations logop not 0;
        # solving a quadratic equation
        if (same(logop, 0)) {
            p0 = 1/(1 + exp(logrr))
        } else {
            p0 = (-(exp(logrr) + 1) * exp(logop) + sqrt(exp(2 * logop) * (exp(logrr) + 
                1)^2 + 4 * exp(logrr + logop) * (1 - exp(logop))))/(2 * exp(logrr) * 
                (1 - exp(logop)))
        }
        p1 = exp(logrr) * p0
    }
    return(c(p0, p1))
} 

#' Calculate risks from log RR and log OP (vectorised)
#'
#' @param logrr log of relative risk
#' 
#' @param logop log of odds product
#'
#' @details The \eqn{log OP} is defined as \eqn{log OP = log[(P(y=1|x=0)/P(y=0|x=0))*(P(y=1|x=1)/P(y=0|x=1))]}. 
#' 
#' @return a matrix \eqn{(P(y=1|x=0),P(y=1|x=1))} with two columns
#' 
#' @examples getProbRR(0,0)
#' 
#' set.seed(0)
#' logrr = rnorm(10,0,1)
#' logop = rnorm(10,0,1)
#' probs = getProbRR(logrr, logop)
#' colnames(probs) = c("P(y=1|x=0)","P(y=1|x=1)")
#' probs
#' 
#' @export
getProbRR = function(logrr, logop = NA) {
    if(is.matrix(logrr) && ncol(logrr) == 2){
        logop = logrr[,2]
        logrr = logrr[,1]
    } else if(is.na(logop) && length(logrr) == 2){
        logop = logrr[2]
        logrr = logrr[1]
    }
    p0 <- ifelse((logop < (-12)) | (logop > 12) | (logrr < (-12)) | (logrr > 12),
                 ## on the boundary South edge: large -ve logrr or (large -ve logop and -ve
                 ## logrr)
                 ifelse ((logrr < (-12)) | ((logop < (-12)) & (logrr < 0)),
                         getPrbAux(logop-logrr),
                         ifelse((logrr > 12) | ((logop < (-12)) & (logrr > 0)),
                                ## West edge: large +ve logrr or (large -ve logop and +ve logrr)
                                0,
                                pmin(exp(-logrr), 1))),
                 ## not on the boundary logop = 0; solving linear equations logop not 0;
                 ## solving a quadratic equation
          ifelse(same(logop, 0),
                 1/(1 + exp(logrr)),
                 (-(exp(logrr) + 1) * exp(logop) + sqrt(exp(2 * logop) * (exp(logrr) + 1)^2 + 4 * exp(logrr + logop) * (1 - exp(logop))))/(2 * exp(logrr) * (1 - exp(logop)))))
    p1 <- ifelse((logop < (-12)) | (logop > 12) | (logrr < (-12)) | (logrr > 12),
                 ## on the boundary South edge: large -ve logrr or (large -ve logop and -ve
                 ## logrr
                 ifelse ((logrr < (-12)) | ((logop < (-12)) & (logrr < 0)),
                         0,
                         ifelse((logrr > 12) | ((logop < (-12)) & (logrr > 0)),
                                ## West edge: large +ve logrr or (large -ve logop and +ve logrr)
                                getPrbAux(logop + logrr),
                                pmin(exp(logrr), 1))),
                 ## not on the boundary logop = 0
                 exp(logrr) * p0)
    cbind(p0,p1)
} 
