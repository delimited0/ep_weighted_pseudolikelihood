#' one factor updating EP spike and slab prior Bayesian variable selection, as in
#' Hernandez-Lobato et al 2013
#' @export
ep_ss1 = function(X, y, v_noise, v_slab, p_incl, 
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
  iter = 0
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
    # Xv = matrix(upsilon, nrow = 1)
    # Xs = X * matrix(v_site2, n, d, byrow = TRUE)  # X \Lambda
    # M = solve(In / v_noise + Xs %*% t(X))  # inv(I * sigma^2 + X \Lambda X^T)
    # 
    # Tv = Xs %*% t(Xv)
    # Cv = t.default(Xv) * matrix(v_site2, d, nrow(Xv)) - t(Xs) %*% (M %*% Tv)
    # nutVnu <- colSums(Cv * t.default(Xv))
    nutVnu = t(upsilon) %*% V %*% upsilon
    
    # logdetV = (sum(log(v_site2)) - determinant(diag(n) + v_noise * X %*% (matrix(v_site2, d, n) * t(X)))$modulus[[1]])
    logdetV = determinant(V)$modulus

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
  
    # Delta_X <- matrix(v_site2, d, n) * t(X)
    # vectMul <- (t(X) %*% m_site2 / v_noise + v_site2^-1 * m_site2)
    # 
    # Inv <- solve(v_noise * In + X %*% (matrix(v_site2, d, n) * X), tol = 1e-100)
    # 
    # # special computation for just the marginals of approximating mean and variance
    # m <- v_site2 * vectMul - Delta_X %*% Inv %*% (t(Delta_X) %*% vectMul)
    # v <- v_site2 - colSums(((Inv) %*% t(Delta_X)) * t(Delta_X))
   
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


#' two factor updating EP spike and slab prior Bayesian variable selection, as in
#' Hernandez-Lobato 2010 and Hernandez-Lobato et al 2015.
#'
#' @param X covariate matrix (n x p)
#' @param y response vector (n x 1)
#' @param v_noise noise variance hyperparameter
#' @param v_slab slab variance hyperparameter
#' @param p_incl inclusion probability hyperparameter
#' @param v_inf infinite site variance stand in
#' @param max_iter maximum number of iterations
#' @param delta convergence parameter change threshold
#' @param k damping multiplier
#' @param damping initial damping value
#' @param woodbury boolean, use woodbury form of update for V or not?
#' @param opt boolean, optimize hyperparameters or not?
#' @export
ep_ss2 = function(X, y, v_noise, v_slab, p_incl, v_inf = 100, max_iter = 200, 
                  delta = 1e-4, k = .99, damping = .5, woodbury = FALSE, 
                  opt = TRUE, opt_method = "Nelder-Mead",
                  lb = c(0, 0), ub = c(Inf, Inf)) {
  
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
  iter = 0
  mlik_value = NA
  
  # define likelihood and gradient functions if optimizing variances
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
  
  if (opt_method == "L-BFGS-B") {
    grad = function(params) {
      v_noise = params[1]
      v_slab = params[2]
      
      cons = 
        plogis(p_site3) * dnorm(0, m_site1, v_site1 + v_slab) +
        plogis(-p_site3) * dnorm(0, m_site1, v_site1)
      
      if (woodbury) {
        dlogalpha_dvnoise = -sum(diag(XV_tX %*% solve(In + XV_tX / v_noise))) / v_noise^2
      }
      else {
        dlogalpha_dvnoise = -sum(diag(V_site2 %*% tXX %*% solve(diag(p) + V_site2 %*% tXX / v_noise))) / v_noise^2    
      }
      dlogs1_dvnoise = -.5 * tmtXy / v_noise^2 - .5 * yty - .5 * n / v_noise
      dlogs2_dvslab = .5 * sum( 
        cons^(-1) * dnorm(0, m_site1, v_site1 + v_slab) * 
          (m_site1^2 - v_site1 - v_slab) / (v_site1 + v_slab^2)
      )
      
      return(c(
        dlogalpha_dvnoise + dlogs1_dvnoise,
        dlogs2_dvslab
      ))
    }
  }
  
  # EP iterations ----
  
  while (!converged & iter < max_iter) {
    
    if (any(v_site1 < 0))
      browser()
    
    # second factor
    p_site2_new = .5*log(v_site1) - .5*log(v_site1 + v_slab) + 
      .5 * m_site1^2 * ( (1 / v_site1) - (1 / (v_site1 + v_slab)) )
    
    a = plogis(p_site2_new + p_site3) * (m_site1 / (v_site1 + v_slab)) +
      plogis(-p_site2_new - p_site3) * (m_site1 / v_site1)
    b = plogis(p_site2_new + p_site3) * ((m_site1^2 - v_site1 - v_slab) / (v_site1 + v_slab)^2) +
      plogis(-p_site2_new - p_site3) * ((m_site1 / v_site1)^2 - (1 / v_site1))
    
    v_site2_new = (1 / (a^2 - b)) - v_site1
    m_site2_new = m_site1 - a * (v_site2_new + v_site1)
    
    # deal with negative variance
    is_neg = v_site2_new < 0
    v_site2_new[is_neg] = v_inf
    
    # overall approximation
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
    p_old = p
    p = p_site2_new + p_site3
    
    # first factor
    v_site1_new = 1 / ( (1 / v) - (1 / v_site2_new) )
    m_site1_new = ( (m / v) - (m_site2_new / v_site2_new) ) * v_site1_new
    
    # damp the updates
    damping = damping * k  
    
    p_site2_damp = damping * p_site2_new + (1-damping) * p_site2
    v_site2_damp = 1 / (damping * (1 / v_site2_new) + (1-damping) * (1 / v_site2))
    m_site2_damp = v_site2_damp * (damping * (m_site2_new / v_site2_new)  + (1-damping) * (m_site2 / v_site2))
    p_site2 = p_site2_damp
    v_site2 = v_site2_damp
    m_site2 = m_site2_damp
    
    v_site1_damp = 1 / (damping * (1 / v_site1_new) + (1 - damping) * (1 / v_site1))
    m_site1_damp = v_site1_damp * (damping * (m_site1_new / v_site1_new) + (1 - damping) * (m_site1 / v_site1))
    v_site1 = v_site1_damp
    m_site1 = m_site1_damp
    
    # check convergence
    if (any(is.na(p))) {
      browser()  
    }
    if (max(abs(p - p_old)) < delta) {
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
    
    # optimize hyperparameters
    if (opt) {
      
      if (opt_method == "Nelder-Mead") {
        hyper_opt = dfoptim::nmkb(par = c(v_noise, v_slab),
                                  fn = mlik,
                                  lower = lb, upper = ub)
      }
      else if (opt_method == "L-BFGS-B") {
        # hyper_opt = optim(par = c(v_noise, v_slab), fn = mlik, gr = grad,
        #                   lower = lb, upper = ub, 
        #                   method = opt_method)  
        hyper_opt = 
          lbfgsb3c::lbfgsb3c(par = c(v_noise, v_slab), fn = mlik, gr = grad,
                             lower = lb, upper = ub)
      }
      else if (opt_method == "BOBYQA") {
        hyper_opt = minqa::bobyqa(par = c(v_noise, v_slab),
                                  fn = mlik,
                                  lower = lb, upper = ub)
      }
      else {
        stop("Invalid optimization method")
      }
      
      
      v_noise = hyper_opt$par[1]
      v_slab = hyper_opt$par[2]
      # mlik_value = -hyper_opt$value
    }
    
    iter = iter + 1
  }
  
  mlik_value = -mlik(c(v_noise, v_slab))
  
  result = list(
    m = m,
    v = v,
    logisp = p,
    p = plogis(p),
    iters = iter,
    v_noise = v_noise,
    v_slab = v_slab,
    llik = mlik_value,
    opt_result = hyper_opt
  )
  
  return(result)
}