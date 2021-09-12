#' @param X covariate matrix (n x p)
#' @param y response vector (n x 1)
#' @param v_noise noise variance hyperparameter
#' @param v_slab slab variance hyperparameter
#' @param p_incl inclusion probability hyperparameter
#' @param v_inf infinite site variance stand in
#' @param max_iter maximum number of iterations
#' @param delta convergence parameter change threshold
#' @param k damping multiplier
#' @param eps initial damping value
#' @param woodbury boolean, use woodbury form of update for V or not?
#' @param opt boolean, optimize hyperparameters or not?
#' @export
ep_wlr = function(X, y, v_noise, v_slab, p_incl, v_inf = 100, max_iter = 200, 
                  delta = 1e-4, k = .99, eps = .5, woodbury = FALSE, opt = TRUE) {
  
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
  tXX = t(X) %*% X
  In = diag(n)
  
  # first updates ----
  
  # third factor
  p_site3 = qlogis(p_incl)
  
  # second factor
  v_site2 = p_incl * v_slab * rep(1, d)
  
  # first factor
  V_site2 = diag(v_site2)
  V_site2_inv = diag(1 / v_site2)
  V = V_site2 - V_site2 %*% t(X) %*%
    solve( v_noise * diag(n) + X %*% V_site2 %*% t(X), X %*% V_site2 )
  v = diag(V)
  m = V %*% ( (1/v_noise) * tXy + V_site2_inv %*% m_site2 )
  v_site1 = 1 / ( (1 / v) - (1 / v_site2) )
  m_site1 = ( (m / v) - (m_site2 / v_site2) ) * v_site1
  
  converged = FALSE
  iter = 1
  mlik_value = NA
  
  # EP iterations ----
  
  while (!converged & iter < max_iter) {
    
    # second factor
    p_site2_new = .5*log(v_site1) - .5*log(v_site1 + v_slab) + 
      .5 * m_site1^2 * ( (1 / v_site1) - (1 / (v_site1 + v_slab)) )
    
    a = plogis(p_site2_new + p_site3) * (m_site1 / (v_site1 + v_slab)) +
      plogis(-p_site2_new - p_site3) * (m_site1 / v_site1)
    b = plogis(p_site2_new + p_site3) * ((m_site1^2 - v_site1 - v_slab) / (v_site1 + v_slab)^2) +
      plogis(-p_site2_new - p_site3) * ((m_site1 / v_site1)^2 - (1 / v_site1))
    
    v_site2_new = (1 / (a^2 - b)) - v_site1
    m_site2_new = m_site1 - a * (v_site2_new + v_site1)
    
    # first factor
    V_site2 = diag(as.vector(v_site2_new))
    V_site2_inv = diag(1 / as.vector(v_site2_new))
    XV_tX = X %*% V_site2 %*% t(X)
    if (woodbury) 
      V = V_site2 - V_site2 %*% t(X) %*% 
      solve( v_noise*In + XV_tX, X %*% V_site2 )
    else {
      V = solve(V_site2_inv + (tXX / v_noise))
    }
    
    v_old = v
    v = diag(V)
    m_old = m
    m = V %*% ( (1/v_noise) * tXy + V_site2_inv %*% m_site2_new)
    p = p_site2_new + p_site3
    
    v_site1_new = 1 / ( (1 / v) - (1 / v_site2_new) )
    m_site1_new = ( (m / v) - (m_site2_new / v_site2_new) ) * v_site1_new
    
    # deal with negative variance
    is_neg <- v_site2_new < 0
    v_site2_new[is_neg] = v_inf
    
    # damp the updates
    eps = eps * k  
    
    p_site2 = eps * p_site2_new + (1 - eps) * p_site2
    v_site2 = 1 / (eps * (1 / v_site2_new) + (1 - eps) * (1 / v_site2))
    m_site2 = v_site2 * (eps * (m_site2_new / v_site2_new)  + (1 - eps) * (m_site2 / v_site2))
    
    v_site1 = 1 / (eps * (1 / v_site1_new) + (1 - eps) * (1 / v_site1))
    m_site1 = v_site1_new * (eps * (m_site1_new / v_site1_new) + (1 - eps) * (m_site1 / v_site1))
    
    # check convergence
    if (all(abs(v - v_old) < delta) & all(abs(m - m_old) < delta)) {
      converged = TRUE
    }
    
    # pre compute marginal likelihood terms
    tmtXy = t(m) %*% tXy
    tmVm = t(m) %*% V_site2_inv %*% m_site2
    tms2Vms2 = t(m_site2) %*% V_site2_inv %*% m_site2
    # sigm_p_site3 = plogis(p_site3, log.p = TRUE)
    # sigm_mp_site3_dnorm = 
    #   plogis(-p_site3, log.p = TRUE) + 
    #   dnorm(0, m_site1, sqrt(v_site1), log = TRUE)
    sum_pmv3 = sum( (m_site1^2 / v_site1) + (m_site2^2 / v_site2) - (m^2 / v) )
    
    mlik = function(params) {
      v_noise = params[1]
      v_slab = params[2]
      
      logs1 = .5*(
        tmVm + (tmtXy / v_noise) - n*log(2*pi * v_noise) - (yty / v_noise) -
          tms2Vms2 -
          determinant(In + (XV_tX / v_noise))$modulus + 
          sum(log1p(v_site2 / v_site1)) + 
          sum_pmv3
      )
      
      logc = log_sum_exp(
        plogis(p_site3, log.p = TRUE) + dnorm(0, m_site1, sqrt(v_site1 + v_slab), log = TRUE),
        plogis(-p_site3, log.p = TRUE) + dnorm(0, m_site1, sqrt(v_site1), log = TRUE)
      )
      
      logs2 = .5*sum(
        2*logc + 
        log1p(v_site1 / v_site2) +
        sum_pmv3 + 
        2*log_sum_exp(
          plogis(p, log.p=TRUE) + plogis(-p_site3, log.p=TRUE),
          plogis(-p, log.p=TRUE) + plogis(p_site3, log.p=TRUE)
        ) -
        2*(plogis(p_site3, log.p=TRUE) + plogis(-p_site3, log.p=TRUE))
      )
      
      value = logs1 + logs2 + .5*d*log(2*pi) +
        .5*sum(log(v)) - .5*sum_pmv3 +
        sum(log_sum_exp(
          plogis(p_site2, log.p = TRUE) + plogis(p_site3, log.p = TRUE),
          plogis(-p_site2, log.p = TRUE) + plogis(-p_site3, log.p = TRUE)
        ))
      
      return(-value)
    }
    
    # optimize hyperparameters
    if (opt) {
      hyper_opt = dfoptim::nmkb(par = c(v_noise, v_slab), 
                                fn = mlik,
                                lower = c(0, 0), upper = c(Inf, Inf))
      
      v_noise = hyper_opt$par[1]
      v_slab = hyper_opt$par[2]
      # mlik_value = -hyper_opt$value
    }
    mlik_value = -mlik(c(v_noise, v_slab))
    
    iter = iter + 1
  }
  
  result <- list(
    m = m,
    v = v,
    logisp = p,
    p = plogis(p),
    iters = iter,
    v_noise = v_noise,
    v_slab = v_slab,
    llik = mlik_value
  )
  
  return(result)
}