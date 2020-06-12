
MLEst = function(param, y, x, va, vb, weights, max.step, thres, alpha.start, 
    beta.start, pa, pb) {
    
    ## starting values for parameter optimization
    if (is.null(alpha.start)) 
        alpha.start = rep(0, pa)
    if (is.null(beta.start)) 
        beta.start = rep(0, pb)
    
    if (param == "OR") {
        fit = stats::glm(y ~ vb - 1 + x * va - va - x, family = "binomial",
                         weights = weights, start = c(beta.start, alpha.start))
        
        point.temp = summary(fit)$coefficients[, 1]
        index = c((pb + 1):(pa + pb), 1:pb)
        point.est = point.temp[index]
        
        cov = stats::vcov(fit)[index, index]
        
        converged = fit$converged
        
    } else {
        
        ### point estimate
        mle = max.likelihood(param, y, x, va, vb, alpha.start, beta.start, 
            weights, max.step, thres, pa, pb)
        point.est = mle$par
        converged = mle$convergence
        # print(point.est)
        alpha.ml = point.est[1:pa]
        beta.ml = point.est[(pa + 1):(pa + pb)]
        
        ### Computing Fisher Information:
        if (param == "RR") 
            cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb, weights)
        if (param == "RD") 
            cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb, weights)
        sd.est = sqrt(diag(cov))
        
    }
    
    name = paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
    sol = WrapResults(point.est, cov, param, name, va, vb, converged)
    return(sol)
    
} 
