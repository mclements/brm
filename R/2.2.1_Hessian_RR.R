Hessian2RR = function(y, x, va, vb, alpha.ml, beta.ml, weights) {
    # calculating the Hessian using the second derivative have to do so
    # because under mis-specification of models Hessian no longer equals the
    # square of the first order derivatives
    
    p0p1 = getProbRR(va %*% alpha.ml, vb %*% beta.ml)
    # p0p1 = cbind(p0, p1): n * 2 matrix
    p0 = p0p1[, 1]
    p1 = p0p1[, 2]
    n = nrow(va)
    pA = p0
    pA[x == 1] = p1[x == 1]
    
    
    ### Building blocks
    
    dpsi0.by.dtheta = -(1 - p0)/(1 - p0 + 1 - p1)
    dpsi0.by.dphi = (1 - p0) * (1 - p1)/(1 - p0 + 1 - p1)
    
    dtheta.by.dalpha = va
    dphi.by.dbeta = vb
    
    dl.by.dpsi0 = (y - pA)/(1 - pA)
    d2l.by.dpsi0.2 = (y - 1) * pA/((1 - pA)^2)
    
    
    
    ###### d2l.by.dalpha.2
    
    d2psi0.by.dtheta.2 = ((p0 - p1) * dpsi0.by.dtheta - (1 - p0) * p1)/((1 - 
        p0 + 1 - p1)^2)
    
    d2l.by.dtheta.2 = d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x)^2 + dl.by.dpsi0 * 
        d2psi0.by.dtheta.2
    
    d2l.by.dalpha.2 = t(dtheta.by.dalpha * d2l.by.dtheta.2 * weights) %*% 
        dtheta.by.dalpha
    
    
    ###### d2l.by.dalpha.dbeta
    
    d2psi0.by.dtheta.dphi = (1 - p0) * (1 - p1) * (p0 - p1)/(1 - p0 + 1 - 
        p1)^3
    
    d2l.by.dtheta.dphi = d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x) * dpsi0.by.dphi + 
        dl.by.dpsi0 * d2psi0.by.dtheta.dphi
    
    d2l.by.dalpha.dbeta = t(dtheta.by.dalpha * d2l.by.dtheta.dphi * weights) %*% 
        dphi.by.dbeta
    d2l.by.dbeta.dalpha = t(d2l.by.dalpha.dbeta)
    # d2l.by.dalpha.dbeta is symmetric itself if (because) va=vb
    
    
    #### d2l.by.dbeta2
    
    d2psi0.by.dphi.2 = (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2)/(1 - p0 + 1 - 
        p1)^2) * dpsi0.by.dphi
    
    d2l.by.dphi.2 = d2l.by.dpsi0.2 * (dpsi0.by.dphi)^2 + dl.by.dpsi0 * d2psi0.by.dphi.2
    
    d2l.by.dbeta.2 = t(dphi.by.dbeta * d2l.by.dphi.2 * weights) %*% dphi.by.dbeta
    
    
    
    hessian = -rbind(cbind(d2l.by.dalpha.2, d2l.by.dalpha.dbeta), cbind(d2l.by.dbeta.dalpha, 
        d2l.by.dbeta.2))
    ### NB Note the extra minus sign here
    
    return(list(hessian = hessian, p0 = p0, p1 = p1, pA = pA, dpsi0.by.dtheta = dpsi0.by.dtheta, 
        dpsi0.by.dphi = dpsi0.by.dphi, dtheta.by.dalpha = dtheta.by.dalpha, 
        dphi.by.dbeta = dphi.by.dbeta, dl.by.dpsi0 = dl.by.dpsi0))
    
} 
