#' vectorized logsumexp function
log_sum_exp = function(u, v) {
  maxuv = pmax(u, v)
  maxuv + log(exp(u - maxuv) + exp(v - maxuv))
}

#' @param X covariate matrix (n x p)
#' @param y response vector (n x 1)
#' @param sigma0 noise standard deviation
#' @param p0 prior inclusion probability
#' @param v_slab slab prior variance
#' @param v_inf infinite site variance stand in
#' @param max_iter maximum number of iterations
#' @param delta convergence parameter change threshold
#' @param k damping multiplier
#' @export
ep_wlr = function(X, y, sigma0, p0, v_slab, v_inf = 100, max_iter = 200, 
                  delta = 1e-4, k = .99) {
  
  d = ncol(X)
  n = length(y)
  
  m = rep(0, d)
  v = rep(Inf, d)
  p = rep(0, d)
  
  m_site2 = rep(0, d)
  p_site2 = rep(0, d)
  s_site = rep(1, 3)
  
  # precomputed values
  tXy = t(X) %*% y
  yty = sum(y^2)
  In = diag(n)
  
  # first updates ----
  
  # third factor
  p_site3 = qlogis(p0)
  
  # second factor
  v_site2 = p0 * v_slab * rep(1, d)
  
  # first factor
  V_site2 = diag(v_site2)
  V_site2_inv = diag(1 / v_site2)
  V = V_site2 - V_site2 %*% t(X) %*%
    solve( sigma0^2 * diag(n) + X %*% V_site2 %*% t(X), X %*% V_site2 )
  # V = 
  v = diag(V)
  m = V %*% ( (1/sigma0^2) * tXy + V_site2_inv %*% m_site2 )
  v_site1 = 1 / ( (1 / v) - (1 / v_site2) )
  m_site1 = ( (m / v) - (m_site2 / v_site2) ) * v_site1
  
  converged = FALSE
  iter = 1
  eps = 1
  
  # EP iterations ----
  
  while (!converged & iter < max_iter) {
    
    # damping factor
    eps = eps * k  
    
    # second factor
    p_site2_undamp = .5*log(v_site1) - .5*log(v_site1 + v_slab) + 
      .5 * m_site1^2 * ( (1 / v_site1) - (1 / (v_site1 + v_slab)) )
    p_site2 = eps * p_site2_undamp + (1 - eps) * p_site2
    
    a = plogis(p_site2 + p_site3) * (m_site1 / (v_site1 + v_slab)) +
      plogis(-p_site2 - p_site3) * (m_site1 / v_site1)
    b = plogis(p_site2 + p_site3) * ((m_site1^2 - v_site1 - v_slab) / (v_site1 + v_slab)^2) +
      plogis(-p_site2 - p_site3) * ((m_site1 / v_site1)^2 - (1 / v_site1))

    v_site2_undamp = (1 / (a^2 - b)) - v_site1
    
    # deal with negative variance
    is_neg <- v_site2_undamp < 0
    v_site2_undamp[is_neg] = v_inf
    
    v_site2 = 1 / (eps * (1 / v_site2_undamp) + (1 - eps) * (1 / v_site2))
    m_site2_undamp = m_site1 - a * (v_site2 + v_site1)
    m_site2 = v_site2 * (eps * (m_site2_undamp / v_site2)  + (1 - eps) * (m_site2 / v_site2))
    
    # first factor
    V_site2 = diag(as.vector(v_site2))
    V_site2_inv = diag(1 / as.vector(v_site2))
    XV_tX = X %*% V_site2 %*% t(X)
    V = V_site2 - V_site2 %*% t(X) %*% 
      solve( sigma0^2*In + XV_tX, X %*% V_site2 )
    
    v_old = v
    v = diag(V)
    m_old = m
    m = V %*% ( (1/sigma0^2) * tXy + V_site2_inv %*% m_site2)
    p = p_site2 + p_site3
    
    v_site1_undamp = 1 / ( (1 / v) - (1 / v_site2) )
    v_site1 = 1 / (eps * (1 / v_site1_undamp) + (1 - eps) * (1 / v_site1))
    m_site1_undamp = ( (m / v) - (m_site2 / v_site2) ) * v_site1
    m_site1 = v_site1 * (eps * (m_site1_undamp / v_site1) + (1 - eps) * (m_site1 / v_site1))
    
    # check convergence
    if (all(abs(v - v_old) < delta) & all(abs(m - m_old) < delta)) {
      converged = TRUE
    }
    
    # optimize hyperparameters
    mv_sites = (m_site1^2 / v_site1) + (m_site2^2 / v_site2) - (m^2 / v)
    log1pvv = log(1 + (v_site2 / v_site1))
    
    f <- function(params) {
      
      sigma0 = params[1]
      v_slab = params[2]
      # 
      # browser()
      
      logs1 = .5 * (
        t(m) %*% (V_site2_inv %*% m_site2 + (tXy / sigma0^2)) -
          n * log(2 * pi * sigma0^2) - (yty / sigma0^2) - t(m_site2) %*% V_site2_inv %*% m_site2 -
          determinant(In + (XV_tX / sigma0^2))$modulus + 
          sum( log1pvv + mv_sites )
      )
      logc = log_sum_exp(
        plogis(p_site3, log.p = TRUE) + dnorm(0, m_site1, v_site1 + v_slab, log = TRUE),
        plogis(-p_site3, log.p = TRUE) + dnorm(0, m_site1, v_site1, log = TRUE)
      )
      # browser()
      logs2 = .5 * sum(
        2*logc + log1pvv + mv_sites + 
          2*log( plogis(p) * plogis(-p_site3) + plogis(-p) * plogis(p_site3) ) -
          2*log( plogis(p_site3) * plogis(-p_site3) )
      )
      value = logs1 + logs2 + .5*d*log(2*pi) + 
        .5 * sum(log(v) - mv_sites) +
        sum( log(plogis(p_site2) * plogis(p_site3) + plogis(-p_site2) * plogis(-p_site3)))
      return(-value)
    }
    
    # hyper_opt = optim(c(sigma0, v_slab), f, method = "Nelder-Mead")
    hyper_opt = dfoptim::nmkb(c(sigma0, v_slab), f,
                              lower = c(0, 0), upper = c(Inf, Inf))

    sigma0 = hyper_opt$par[1]
    v_slab = hyper_opt$par[2]
        
    iter = iter + 1
  }
  
  result <- list(
    m = m,
    v = v,
    p = p,
    iters = iter,
    sigma0 = sigma0,
    v_slab = v_slab
  )
  
  return(result)
}