#' @export
ep_ss = function(X, y, v_noise, v_slab, p_incl, 
                 delta = 1e-4, max_iter = 200,
                 damping = .90, k = .99, 
                 v_inf = 100, 
                 woodbury = FALSE, opt = FALSE) {
  
  # initialization
  
  d = ncol(X)
  n = length(y)
  
  m = rep(0, d)
  v = rep(Inf, d)
  p = rep(0, d)
  
  m_site2 = rep(0, d)
  v_site2 = p_incl * v_slab * rep(1, d)
  p_site2 = rep(0, d)
  s_site = rep(1, 3)
  
  Xty = t(X) %*% y
  yty = sum(y^2)
  XtX = t(X) %*% X
  In = diag(n)
  
  # first factor
  Prec_site1 = XtX / v_noise
  mPrec_site1 = Xty / v_noise
  
  # third factor
  p_site3 = qlogis(p_incl)
  
  converged = FALSE
  iter = 1
  mlik_value = NA
  
  # initial approximate parameter
  
  Lambda = diag(v_site2)
  Lambda_inv = diag(1 / v_site2)
  
  eta = mPrec_site1 + Lambda_inv %*% m_site2
  if (woodbury) {
    V = Lambda - Lambda %*% t(X) %*% solve((In / v_noise) + X %*% Lambda %*% t(X)) %*% X %*% Lambda
  }
  else {
    V = solve(Lambda_inv + Prec_site1)
  }
  m = as.vector(V %*% eta)
  v = diag(V)
  
  # negative marginal likelihood function
  neg_mlik = function(params) {
    v_noise = params[1]
    v_slab = params[2]

    logk = log_sum_exp(
      plogis(p_cavity, log.p=TRUE) + dnorm(0, m_cavity, sqrt(v_slab + v_cavity), log=TRUE),
      plogis(-p_cavity, log.p=TRUE) + dnorm(0, m_cavity, sqrt(v_cavity), log=TRUE)
    )

    logs1 = -.5 * n * log(2*pi * v_noise) - .5 * yty / v_noise
    logs2 = sum(
      logk -
      log_sum_exp(
        plogis(p_site2, log.p=TRUE) + plogis(p_cavity, log.p=TRUE),
        plogis(-p_site2, log.p=TRUE) + plogis(-p_cavity, log.p=TRUE)) -
      .5 * log(2*pi * v_site2) -
      dnorm(0, m_cavity - m_site2, sqrt(v_site2 + v_cavity), log=TRUE)
    )

    upsilon = mPrec_site1 + Lambda_inv %*% m_site2
    Xv = matrix(upsilon, nrow = 1)
    Xs = X * matrix(v_site2, n, d, byrow = TRUE)  # X \Lambda
    M = solve(In / v_noise + Xs %*% t(X))  # inv(I * sigma^2 + X \Lambda X^T)
    
    Tv = Xs %*% t(Xv)
    Cv = t.default(Xv) * matrix(v_site2, d, nrow(Xv)) - t(Xs) %*% (M %*% Tv)
    nutVnu <- colSums(Cv * t.default(Xv))
    
    logdetV = (sum(log(v_site2)) - determinant(diag(n) + v_noise * X %*% (matrix(v_site2, d, n) * t(X)))$modulus[[1]])

    value =
      logs1 + logs2 +
      sum(log_sum_exp(
        plogis(p_site3, log.p=TRUE) + sum(plogis(p_site2, log.p=TRUE)),
        plogis(-p_site3, log.p=TRUE) + sum(plogis(-p_site2, log.p=TRUE))
      )) -
      .5 * m_site2 %*% Lambda_inv %*% m_site2 + 
      .5 * d*log(2*pi) +
      .5 * logdetV +
      .5 * nutVnu

    return(-value)
  }

  # second factor updates ----
  while (!converged & iter <= max_iter) {
    
    # cavity parameters
    v_cavity = (v^(-1) - v_site2^(-1))^(-1)
    m_cavity = v_cavity * (v^(-1) * m - v_site2^(-1) * m_site2)
    p_cavity = p - p_site2
    
    # second factor parameters
    p_site2_new = 
      dnorm(0, m_cavity, sqrt(v_cavity + v_slab), log=TRUE) -
      dnorm(0, m_cavity, sqrt(v_cavity), log=TRUE)
    
    a = 
      plogis(p_site2_new + p_cavity) * (m_cavity / (v_cavity + v_slab)) + 
      plogis(-p_site2_new - p_cavity) * (m_cavity / v_cavity)
    b = 
      plogis(p_site2_new + p_cavity) * (m_cavity^2 - v_cavity - v_slab) / (v_cavity + v_slab)^2 +
      plogis(-p_site2_new - p_cavity) * (m_cavity^2 - v_cavity) / v_cavity^2
    
    v_site2_new = (a^2 - b)^(-1) - v_cavity
    m_site2_new = m_cavity + a / (a^2 - b)
    
    # handle negative variances
    update_dims = !(v_cavity < 0)
    update_dims = 1:d
    v_site2_negative = v_site2_new < 0
    v_site2_new[v_site2_negative] = v_inf
    
    # damped update
    p_site2_old = p_site2
    v_site2_old = v_site2
    m_site2_old = m_site2
    
    p_site2[update_dims] = 
      damping * p_site2_new[update_dims] + (1-damping) * p_site2_old[update_dims]
    v_site2[update_dims] = 
      ( damping * v_site2_new[update_dims]^(-1) + (1-damping) * v_site2_old[update_dims]^(-1) )^(-1)
    m_site2[update_dims] = 
      v_site2 * (damping * m_site2_new / v_site2_new + (1-damping) * m_site2_old / v_site2_old)
    
    p_old = p
    p = p_site2 + p_site3
    
    Lambda = diag(v_site2)
    Lambda_inv = diag(1 / v_site2)
    
    eta = mPrec_site1 + Lambda_inv %*% m_site2
    if (woodbury) {
      V = Lambda - Lambda %*% t(X) %*% solve((In / v_noise) + X %*% Lambda %*% t(X)) %*% X %*% Lambda
    }
    else {
      V = solve(Lambda_inv + Prec_site1)
    }
    m = as.vector(V %*% eta)
    v = diag(V)
   
    # check convergence
    if (max(abs(p - p_old)) < delta)
      converged = TRUE
    
    # optimize hyperparameters
    if (opt) {
    }
    
    iter = iter + 1 
    damping = damping * k
  }
  
  if (!opt) {
    mlik_value = -neg_mlik(c(v_noise, v_slab))
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