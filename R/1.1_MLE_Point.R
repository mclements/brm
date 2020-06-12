
max.likelihood = function(param, y, x, va, vb, alpha.start, beta.start, weights, 
                          max.step, thres, pa, pb) {
    
    startpars = c(alpha.start, beta.start)
    
    getProb = if (param == "RR") getProbRR else getProbRD
    
    ## negative log likelihood function
    neg.log.likelihood = function(pars) {
        alpha = pars[1:pa]
        beta = pars[(pa + 1):(pa + pb)]
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0 = p0p1[, 1];   p1 = p0p1[, 2]
        
        return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weights[x == 0] + 
                        (y[x == 0]) * log(p0[x == 0]) * weights[x == 0]) - sum((1 - y[x == 
                                                                                          1]) * log(1 - p1[x == 1]) * weights[x == 1] + (y[x == 1]) * log(p1[x == 
                                                                                                                                                                 1]) * weights[x == 1]))
    }
    
    neg.log.likelihood.alpha = function(alpha){
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0    = p0p1[,1];  p1 = p0p1[,2]
        
        return(-sum((1-y[x==0])*log(1-p0[x==0])*weights[x==0] +
                        (y[x==0])*log(p0[x==0])*weights[x==0]) -
                   sum((1-y[x==1])*log(1-p1[x==1])*weights[x==1] +
                           (y[x==1])*log(p1[x==1])*weights[x==1]))  
    }
    
    neg.log.likelihood.beta = function(beta){
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0    = p0p1[,1];  p1 = p0p1[,2]
        
        return(-sum((1-y[x==0])*log(1-p0[x==0])*weights[x==0] +
                        (y[x==0])*log(p0[x==0])*weights[x==0]) -
                   sum((1-y[x==1])*log(1-p1[x==1])*weights[x==1] +
                           (y[x==1])*log(p1[x==1])*weights[x==1]))  
    }
    
    
    ## Optimization 
    
    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    alpha = alpha.start; beta = beta.start
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
        step = step + 1
        opt1 = stats::optim(alpha,neg.log.likelihood.alpha,control=list(maxit=max(100,max.step/10)))
        diff1 = Diff(opt1$par,alpha)
        alpha = opt1$par
        opt2 = stats::optim(beta,neg.log.likelihood.beta,control=list(maxit=max(100,max.step/10)))
        diff  = max(diff1,Diff(opt2$par,beta))
        beta = opt2$par
    }
    
    opt = list(par = c(alpha,beta), convergence = (step < max.step), 
               value = neg.log.likelihood(c(alpha,beta)), step = step)
    
    return(opt)
}

