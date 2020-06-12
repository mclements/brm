

### variance calculation

var.mle.rr = function(x, alpha.ml, beta.ml, va, vb, weights) {
    
    p0p1 = getProbRR(va %*% alpha.ml, vb %*% beta.ml)
    n = dim(va)[1]
    pA = rep(NA, n)     # P(Y=1|A,V); here A = X
    pA[x == 0] = p0p1[x == 0, 1]
    pA[x == 1] = p0p1[x == 1, 2]
    
    expect.dl.by.dpsi0.squared = (pA)/(1 - pA)
    dpsi0.by.dphi = (1 - p0p1[, 1]) * (1 - p0p1[, 2])/((1 - p0p1[, 1]) + (1 - 
        p0p1[, 2]))
    dpsi0.by.dtheta = -(1 - p0p1[, 1])/((1 - p0p1[, 1]) + (1 - p0p1[, 2]))
    tmp = cbind((dpsi0.by.dtheta + x) * va, dpsi0.by.dphi * vb)
    ## since dtheta.by.dalpha = va, and dphi.by.dbeta = vb
    fisher.info = (t(expect.dl.by.dpsi0.squared * weights * tmp) %*% tmp)
    return(solve(fisher.info))
}




### variance calculation

var.mle.rd = function(x, alpha.ml, beta.ml, va, vb, weights) {
    
    p0p1 = getProbRD(va %*% alpha.ml, vb %*% beta.ml)
    # p0p1 = cbind(p0, p1): n * 2 matrix
    p0 = p0p1[, 1]
    p1 = p0p1[, 2]
    n = nrow(va)
    pA = p0             # P(Y=1|A,V); here A = X
    pA[x == 1] = p1[x == 1]
    s0 = p0 * (1 - p0)
    s1 = p1 * (1 - p1)
    sA = pA * (1 - pA)
    
    rho = as.vector(tanh(va %*% alpha.ml))  #estimated risk differences
    
    expect.dl.by.dpA.squared = 1/sA
    dp0.by.dphi = s0 * s1/(s0 + s1)
    dp0.by.drho = -s0/(s0 + s1)
    drho.by.dalpha = (1 - rho^2) * va
    dphi.by.dbeta = vb
    
    tmp = cbind((dp0.by.drho + x) * drho.by.dalpha, dp0.by.dphi * dphi.by.dbeta)
    fisher.info = (t(expect.dl.by.dpA.squared * weights * tmp) %*% tmp)
    return(solve(fisher.info))
} 
