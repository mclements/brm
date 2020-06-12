
DREst = function(param, y, x, va, vb, vc, alpha.ml, beta.ml, gamma, optimal, 
    weights, max.step, thres, alpha.start, beta.cov, gamma.cov, message) {
    
    dr.est = dr.estimate.noiterate(param, y, x, va, vb, vc, alpha.ml, beta.ml, 
        gamma, optimal, weights, max.step, thres, alpha.start, message)
    point.est = dr.est$par
    converged = dr.est$convergence
    
    if (param == "RR") 
        alpha.cov = var.rr.dr(y, x, va, vb, vc, point.est, alpha.ml, beta.ml, 
                              gamma, optimal, weights)
    if (param == "RD") 
        alpha.cov = var.rd.dr(y, x, va, vb, vc, point.est, alpha.ml, beta.ml,
                              gamma, optimal, weights)
    
    pa = dim(va)[2]
    pb = dim(vb)[2]
    pc = dim(vc)[2]
    name = paste(c(rep("alpha", pa), rep("beta", pb), rep("gamma", pc)),
                 c(1:pa, 1:pb, 1:pc))
    point.est = c(point.est, beta.ml, gamma)
    cov = matrix(NA,pa+pb+pc, pa+pb+pc)
    cov[1:pa,1:pa] = alpha.cov
    cov[(pa+1):(pa+pb),(pa+1):(pa+pb)] = beta.cov
    cov[(pa+pb+1):(pa+pb+pc),(pa+pb+1):(pa+pb+pc)] = gamma.cov
                 
    sol = WrapResults(point.est, cov, param, name, va, vb, converged)
    return(sol)
    
} 
