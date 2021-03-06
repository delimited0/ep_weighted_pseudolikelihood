#' @param X covariate matrix (n x p)
#' @param y response vector (n x 1)
#' @param sigma0 noise standard deviation
#' @param p0 prior inclusion probability
#' @param v_slab slab prior variance
#' @param v_inf infinite site variance stand in
#' @param max_iter maximum number of iterations
#' @param delta convergence parameter change threshold
#' @param k damping multiplier
#' @param woodbury boolean, use woodbury form of update for V or not?
#' @param opt boolean, optimize hyperparameters or not?
#' @export
ep_wlr_nmp = function(X, y, sigma0, p0, v_slab, v_inf = 100, max_iter = 200, 
                  delta = 1e-4, k = .99, woodbury = FALSE, opt = TRUE,
                  p0_limit = .2) {
  
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
  p_site3 = qlogis(p0)
  
  # second factor
  v_site2 = rep(p0 * v_slab, d)
  
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
  mlik_value = NA
  
  # EP iterations ----
  
  while (!converged & iter < max_iter) {
    
    # damping factor
    eps = eps * k  
    
    # first factor 
    p_site3 = qlogis(p0)
    
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
    if (woodbury) 
      V = V_site2 - V_site2 %*% t(X) %*% 
      solve( sigma0^2*In + XV_tX, X %*% V_site2 )
    else {
      V = solve(V_site2_inv + (tXX / sigma0^2))
    }
    
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
    # mv_sites = (m_site1^2 / v_site1) + (m_site2^2 / v_site2) - (m^2 / v)
    # log1pvv = log(1 + (v_site2 / v_site1))
    tmtXy = t(m) %*% tXy
    sigm_p_site3 = plogis(p_site3, log.p = TRUE)
    sigm_mp_site3_dnorm = plogis(-p_site3, log.p = TRUE) + dnorm(0, m_site1, v_site1, log = TRUE)
    
    if (opt) {
      mlik = function(params) {
        sigma0 = params[1]
        v_slab = params[2]
        p0 = params[3]
        
        logs1 = .5*(
          tmtXy / sigma0^2 - n * log(sigma0^2) - (yty / sigma0^2) -
            determinant(In + (XV_tX / sigma0^2))$modulus
        )
        
        logc = log_sum_exp(
          # sigm_p_site3 + dnorm(0, m_site1, v_site1 + v_slab, log = TRUE),
          # sigm_mp_site3_dnorm
          log(p0) + dnorm(0, m_site1, v_site1 + v_slab, log = TRUE),
          log(1 - p0) + dnorm(0, m_site1, v_site1, log = TRUE)
        )
        
        logps = 
          log_sum_exp(
            plogis(p_site2 + qlogis(p0), log.p = TRUE) + log(1-p0),
            plogis(-p_site2 - qlogis(p0), log.p = TRUE) + log(p0)
            # plogis(p, log.p=TRUE) + log(1-p0),
            # plogis(-p, log.p=TRUE) + log(p0)
          ) -
          (log(p0) + log(1-p0))
        logs2 = sum(logc + logps) 
        
        logep = sum(log_sum_exp(
          plogis(p_site2, log.p = TRUE) + log(p0),
          plogis(-p_site2, log.p = TRUE) + log(1-p0)
        ))
        
        value = logs1 + logs2 + logep
        
        return(-value)
      }
      hyper_opt = dfoptim::nmkb(par = c(sigma0, v_slab, p0), 
                                fn = mlik,
                                lower = c(0, 0, 0), upper = c(Inf, Inf, p0_limit))
      
      # fn = mlik_obj,
      # n = n, tmtXy = tmtXy, yty = yty, In = In,
      # XV_tX = XV_tX,
      # sigm_p_site3 = sigm_p_site3,
      # sigm_mp_site3_dnorm = sigm_mp_site3_dnorm,
      # m_site1 = m_site1, v_site1 = v_site1)
      
      sigma0 = hyper_opt$par[1]
      v_slab = hyper_opt$par[2]
      p0 = hyper_opt$par[3]
      mlik_value = -hyper_opt$value
    }
    
    iter = iter + 1
  }
  
  result <- list(
    m = m,
    v = v,
    p = p,
    iters = iter,
    sigma0 = sigma0,
    v_slab = v_slab,
    p0 = p0,
    llik = mlik_value
  )
  
  return(result)
}