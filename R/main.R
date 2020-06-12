#' Fitting Binary Regression Models
#'
#' @description  \code{brm} is used to estimate the association between two binary variables, and how that varies as a function of other covariates.
#'
#' @param param The measure of association. Can take value 'RD' (risk difference),
#'  'RR' (relative risk) or 'OR' (odds ratio)
#'  
#' @param y The response vector. Should only take values 0 or 1.
#' 
#' @param x The exposure vector. Should only take values 0 or 1.
#' 
#' @param va The covariates matrix for the target model. It can be specified via an object of class "\code{formula}" or a matrix. In the latter case, no intercept terms will be added to the covariates matrix.
#' 
#' @param vb The covariates matrix for the nuisance model. It can be specified via an object of class "\code{formula}" or a matrix. In the latter case, no intercept terms will be added to the covariates matrix. (If not specified, defaults to va.)
#' 
#' @param vc The covariates matrix for the probability of exposure, often called the propensity score. It can be specified via an object of class "\code{formula}" or a matrix. In the latter case, no intercept terms will be added to the covariates matrix. By default we fit a logistic regression model for the propensity score. (If not specified, defaults to va.)
#' 
#' @param weights An optional vector of 'prior weights' to be used in the fitting process. Should be NULL or a numeric vector.
#' 
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' 
#' @param est.method The method to be used in fitting the model. Can be 'MLE' (maximum likelihood estimation, the default) or 'DR' (doubly robust estimation).
#' 
#' @param optimal   Use the optimal weighting function for the doubly robust estimator? Ignored if the estimation method is 'MLE'. The default is TRUE.
#' 
#' @param max.step  The maximal number of iterations to be passed into the \code{\link[stats]{optim}} function. The default is 1000.
#' 
#' @param thres Threshold for judging convergence. The default is 1e-6.
#' 
#' @param alpha.start   Starting values for the parameters in the target model.
#' 
#' @param beta.start    Starting values for the parameters in the nuisance model.
#'
#' @param message   Show optimization details? Ignored if the estimation method is 'MLE'. The default is FALSE.
#' 
#' @details  \code{brm} contains two parts: the target model for the dependence measure (RR, RD or OR) and the nuisance model; the latter is required for maximum likelihood estimation. 
#'If \code{param="RR"} then the target model is \eqn{log RR(va) = \alpha*va}. 
#'If \code{param="RD"} then the target model is \eqn{atanh RD(va) = \alpha*va}. 
#'If \code{param="OR"} then the target model is \eqn{log OR(va) = \alpha*va}. 
#'For RR and RD, the nuisance model is for the log Odds Product: \eqn{log OP(vb) = \beta*vb}. 
#'For OR, the nuisance model is for the baseline risk: \eqn{logit(P(y=1|x=0,vb)) = \beta*vb.}
#' In each case the nuisance model is variation independent of the target model, which  ensures that the predicted probabilities lie in \eqn{[0,1]}. 
#' See Richardson et al. (2016+) for more details.
#'
#' If \code{est.method="DR"} then given a correctly specified logistic regression model for the propensity score \eqn{logit(P(x=1|vc)) = \gamma*vc}, estimation of the RR or RD is consistent, even if the log Odds Product model is misspecified. This estimation method is not available for the OR. See Tchetgen Tchetgen et al. (2014) for more details.
#' 
#' When estimating RR and RD, \code{est.method="DR"} is recommended unless it is known that the log Odds Product model is correctly specified. Optimal weights (\code{optimal=TRUE}) are also recommended to increase efficiency.
#' 
#' For the doubly robust estimation method, MLE is used to obtain preliminary estimates of \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}. The estimate of \eqn{\alpha} is then updated by solving a doubly-robust estimating equation. (The estimate for \eqn{\beta} remains the MLE.)
#' 
#' @return A list consisting of
#'  \item{param}{the measure of association.} 
#'
#' \item{point.est}{ the point estimates.}
#' 
#' \item{se.est}{the standard error estimates.}
#' 
#' \item{cov}{estimate of the covariance matrix for the estimates.}
#' 
#' \item{conf.lower}{ the lower limit of the 95\% (marginal) confidence interval.}
#' 
#' \item{conf.upper}{ the upper limit of the 95\% (marginal) confidence interval.}
#' 
#' \item{p.value}{the two sided p-value for testing zero coefficients.}
#' 
#' \item{coefficients}{ the matrix summarizing key information: point estimate, 95\% confidence interval and p-value.}
#' 
#' \item{param.est}{the fitted RR/RD/OR.}
#' 
#' \item{va}{   the matrix of covariates for the target model.}
#' 
#' \item{vb}{   the matrix of covariates for the nuisance model.}
#' 
#' \item{converged}{  Did the maximization process converge? }
#' 
#' @author Linbo Wang <linbo.wang@utoronto.ca>, Mark Clements <mark.clements@ki.se>, Thomas Richardson <thomasr@uw.edu>
#' 
#' @references Thomas S. Richardson, James M. Robins and Linbo Wang. "On Modeling and Estimation for the Relative Risk and Risk Difference." Journal of the American Statistical Association: Theory and Methods (2017).
#' 
#' Eric J. Tchetgen Tchetgen, James M. Robins and Andrea Rotnitzky. "On doubly robust estimation in a semiparametric odds ratio model." Biometrika 97.1 (2010): 171-180.
#' 
#' @seealso \code{getProbScalarRD}, \code{getProbRD} (vectorised), \code{getProbScalarRR} and \code{getProbRR} (vectorised) for functions calculating risks P(y=1|x=1) and P(y=1|x=0) from (atanh RD, log OP) or (log RR, log OP);
#' 
#' \code{predict.blm} for obtaining fitted probabilities from \code{brm} fits.
#' 
#' @examples 
#' set.seed(0)
#' n           = 100
#' alpha.true  = c(0,-1)
#' beta.true   = c(-0.5,1)
#' gamma.true  = c(0.1,-0.5)
#' params.true = list(alpha.true=alpha.true, beta.true=beta.true, 
#'      gamma.true=gamma.true)
#' v.1         = rep(1,n)       # intercept term
#' v.2         = runif(n,-2,2) 
#' v           = cbind(v.1,v.2)
#' pscore.true = exp(v %*% gamma.true) / (1+exp(v %*% gamma.true))
#' p0p1.true   = getProbRR(v %*% alpha.true,v %*% beta.true)
#' x           = rbinom(n, 1, pscore.true)  
#' pA.true       = p0p1.true[,1]
#' pA.true[x==1] = p0p1.true[x==1,2]
#' y = rbinom(n, 1, pA.true)
#' 
#' fit.mle = brm(y,x,v,v,'RR','MLE',v,TRUE)
#' fit.drw = brm(y,x,v,v,'RR','DR',v,TRUE)
#' fit.dru = brm(y,x,v,v,'RR','DR',v,FALSE)
#'
#' fit.mle2 = brm(y,x,~v.2, ~v.2, 'RR','MLE', ~v.2,TRUE) # same as fit.mle
#'
#' @export


brm = function(y, x, va, vb = NULL, param, est.method = "MLE", vc = NULL, 
    optimal = TRUE, weights = NULL, subset = NULL, max.step = NULL, thres = 1e-8, 
    alpha.start = NULL, beta.start = NULL, message = FALSE) {
    
    # default param = 'RR'; est.method = 'MLE'; va = v; vb = v; vc = v;
    # weights = NULL; subset = NULL; optimal = TRUE; max.step = NULL; 
    # thres = 1e-06; alpha.start = NULL; beta.start = NULL
    
    if (is.null(vb)) 
        vb = va
    if (is.null(vc)) 
        vc = va
    
    if(class(va)[1] == "formula") va = stats::model.matrix(va)
    if(class(vb)[1] == "formula") vb = stats::model.matrix(vb)
    if(class(vc)[1] == "formula") vc = stats::model.matrix(vc)
    
    if(is.vector(va)) va = as.matrix(va, ncol = 1)
    if(is.vector(vb)) vb = as.matrix(vb, ncol = 1)
    if(is.vector(vc)) vc = as.matrix(vc, ncol = 1)
    
    if (is.null(weights)) 
        weights = rep(1, length(y))
    if (is.null(subset))
        subset = 1:length(y)
    
    ValidCheck(param, y, x, va, vb, vc, weights, subset, est.method, optimal,
               max.step, thres, alpha.start, beta.start)
    
    data = cbind(y,x,va,vb,vc,weights)[subset,]
    subset = subset[rowSums(is.na(data)) == 0]
    y = y[subset];  x = x[subset];  va = va[subset,,drop=FALSE];    vb = vb[subset,,drop=FALSE]
    vc = vc[subset,,drop=FALSE];   weights = weights[subset]
    
    pa = dim(va)[2]
    pb = dim(vb)[2]
    if (is.null(max.step)) max.step = min(pa * 20, 1000)
    
    if (est.method == "MLE"){
        sol = MLEst(param, y, x, va, vb, weights, max.step, thres, alpha.start, 
                    beta.start, pa, pb)
    }
    if (est.method == "DR") {
        if (param == "OR") {
            cat("No doubly robust estimation methods for OR (with propensity score models) are available. Please refer to Tchetgen Tchetgen et al. (2010) for an alternative doubly robust estimation method. \n")
            return()
        }
        if (is.null(alpha.start) | is.null(beta.start)){
            sol = MLEst(param, y, x, va, vb, weights, max.step, thres, 
                        alpha.start, beta.start, pa, pb)
            alpha.ml = sol$point.est[1:pa]
            beta.ml = sol$point.est[(pa+1):(pa+pb)]
            beta.cov = sol$cov[(pa+1):(pa+pb),(pa+1):(pa+pb)]
            alpha.start = alpha.ml
        }else{
            alpha.ml = alpha.start; beta.ml = beta.start
            beta.cov = matrix(NA,pb,pb)
        } 
        
        gamma.fit = stats::glm(x~vc-1, weight = weights, family="binomial")
        gamma = gamma.fit$coefficients
        gamma.cov = summary(gamma.fit)$cov.unscaled
        sol = DREst(param, y, x, va, vb, vc, alpha.ml, beta.ml, gamma, optimal, 
            weights, max.step, thres, alpha.start, beta.cov, gamma.cov, message)
    }
    
    return(sol)
    
} 
