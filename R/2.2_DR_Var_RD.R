
## Sandwich estimator for variance of RD

var.rd.dr = function(y, x, va, vb, vc, alpha.dr, alpha.ml, beta.ml, gamma, 
    optimal, weights) {
    ######################################## 
    
    pscore = as.vector(expit(vc %*% gamma))
    n = length(pscore)
    
    ### 1. - E[dS/d(alpha.ml,beta.ml)] ############################## Computing
    ### the Hessian:
    
    Hrd = Hessian2RD(y, x, va, vb, alpha.ml, beta.ml, weights)
    hessian = Hrd$hessian
    p0  = Hrd$p0; p1 = Hrd$p1; pA = Hrd$pA
    s0  = Hrd$s0; s1 = Hrd$s1; sA = Hrd$sA
    rho = Hrd$rho
    dl.by.dpA      = Hrd$dl.by.dpA
    dp0.by.dphi    = Hrd$dp0.by.dphi
    dp0.by.drho    = Hrd$dp0.by.drho
    drho.by.dalpha = Hrd$drho.by.dalpha
    dphi.by.dbeta  = Hrd$dphi.by.dbeta
    dpA.by.drho    = Hrd$dpA.by.drho
    dpA.by.dalpha  = Hrd$dpA.by.dalpha
    dpA.by.dphi    = Hrd$dpA.by.dphi
    dpA.by.dbeta   = Hrd$dpA.by.dbeta
    
    
    ############# extra building blocks ##########################
    
    H.alpha = y - x * as.vector(tanh(va %*% alpha.dr))
    
    ############# Calculation of optimal vector (used in several places below) ##
    
    if (optimal == TRUE) {
        wt = (1 - rho^2)/(pscore * s0 + (1 - pscore) * s1)
    } else {
        wt = rep(1, n)
    }
    
    
    ### 2. -E[dU.by.dalphaml.betaml] ####################################
    
    dU.by.dp0 = -va * wt * (x - pscore)  # n by 2
    dp0.by.dalpha.beta = cbind(drho.by.dalpha * dp0.by.drho, dphi.by.dbeta * 
        dp0.by.dphi)  # n by 4
    
    dU.by.dwt = va * (x - pscore) * (H.alpha - p0)  # n by 2
    
    if (optimal == TRUE) {
        esA = pscore * s0 + (1 - pscore) * s1  # E[s_{1-A}]...
        dwt.by.drho = (-2 * rho * esA - (1 - rho^2) * (1 - pscore) * (1 - 
            2 * p1))/esA^2
        dwt.by.dp0 = -(1 - rho^2) * (2 * pscore * rho + 1 - 2 * p1)/esA^2
        
        dwt.by.dalpha = drho.by.dalpha * (dwt.by.drho + dwt.by.dp0 * dp0.by.drho)
        dwt.by.dbeta = dphi.by.dbeta * dwt.by.dp0 * dp0.by.dphi
        
        dwt.by.dalpha.beta = cbind(dwt.by.dalpha, dwt.by.dbeta)  # n by 4    
    } else {
        dwt.by.dalpha.beta = matrix(0, n, ncol(va) + ncol(vb))
    }
    
    dU.by.dalpha.ml.beta.ml = t(dU.by.dp0 * weights) %*% (dp0.by.dalpha.beta) + 
        t(dU.by.dwt * weights) %*% dwt.by.dalpha.beta
    
    
    ### 3. tau = -E[dU/dalpha.dr] ######################################## (This
    ### is the bread of the sandwich estimate)
    
    dU.by.dH = va * wt * (x - pscore)  # n by 2
    rho.dr = as.vector(tanh(va %*% alpha.dr))
    dH.by.dalpha.dr = -va * x * (1 - rho.dr^2)
    
    tau = -t(dU.by.dH * weights) %*% dH.by.dalpha.dr/sum(weights)  # 2 by 2
    
    
    ### 4. E[d(prop score score equation)/dgamma]
    
    dpscore.by.dgamma = vc * pscore * (1 - pscore)  # n by 2   
    part4 = -t(vc * weights) %*% dpscore.by.dgamma  # 2 by 2
    
    
    ### 5. E[dU/dgamma]
    
    dU.by.dpscore = -va * wt * (H.alpha - p0)  # n by 2
    
    if (optimal == TRUE) {
        dwt.by.dpscore = -(1 - rho^2) * (s0 - s1)/esA^2
        dwt.by.dgamma = dpscore.by.dgamma * dwt.by.dpscore  # n by 2
    } else {
        dwt.by.dgamma = matrix(0, n, ncol(vc))
    }
    
    dU.by.dgamma = t(dU.by.dpscore * weights) %*% dpscore.by.dgamma + t(dU.by.dwt * 
        weights) %*% dwt.by.dgamma  # 2 by 2
    
    
    
    ############################################################################# Assembling semi-parametric variance matrix
    
    U = va * wt * (x - pscore) * (H.alpha - p0)  # n by 2
    
    S = cbind(dpA.by.dalpha, dpA.by.dbeta) * dl.by.dpA
    
    pscore.score = vc * (x - pscore)
    
    Utilde = U - t(dU.by.dalpha.ml.beta.ml %*% (-solve(hessian)) %*% t(S)) - 
        t(dU.by.dgamma %*% (solve(part4)) %*% t(pscore.score))  # n by 2
    USigma = t(Utilde * weights) %*% Utilde/sum(weights)
    
    
    
    ################################### Asymptotic var matrix for alpha.dr
    
    alpha.dr.variance = solve(tau) %*% USigma %*% solve(tau)/sum(weights)
    return(alpha.dr.variance)
    
} 
