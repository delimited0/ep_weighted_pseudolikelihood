
grad = function(u, phi, alpha) {
  
  sigmoid_phi = plogis(phi)
  umask = u > (1 - sigmoid_phi)
  dmask = u < sigmoid_phi
  
  
}

log_alpha_div = function(z, phi, X, y, v_noise, v_slab, p_incl) {
  p = ncol(X)
  n = nrow(X)
  
  X_select = X[, z]
  
  cov_mat = v_slab * tcrossprod(X_select) + v_noise * diag(n)
  
  # browser()
  
  llik = 
    (1-alpha) * (
      mvtnorm::dmvnorm(y, rep(0, n), cov_mat, log = TRUE) +
        sum(dbinom(z, 1, p_incl, log = TRUE)) -
        sum(dbinom(z, 1, plogis(phi), log = TRUE))
    )
  
  return(llik)
}

alpha_div = function(z, phi, X, y, v_noise, v_slab, p_incl) {
  p = ncol(X)
  n = nrow(X)
  
  X_select = X[, z]
  
  cov_mat = v_slab * tcrossprod(X_select) + v_noise * diag(n)
  
  # browser()
  
  log_div = 
    (1-alpha) * (
      mvtnorm::dmvnorm(y, rep(0, n), cov_mat, log = TRUE) +
        sum(dbinom(z, 1, p_incl, log = TRUE)) -
        sum(dbinom(z, 1, plogis(phi), log = TRUE))
    )
  
  return(exp(log_div))
}

#' @export
alpha_ss = function(X, y, alpha, v_noise, v_slab, p_incl, 
                    n_div = 30, n_u2g = 50, min_logval = -100,
                    lr = .1, delta = 1e-4, max_iter = 200) {
  
  p = ncol(X)
  n = nrow(X)
  phi = rep(0, p)
  
  converged = FALSE
  iter = 1
  alpha_div_history = rep(NA, max_iter)
  
  # f = function(z) {
  #   cov_mat = X %*% ( (v_slab*diag(p)) * tcrossprod(z) ) %*% t(X) + v_noise*diag(n)
  #   
  #   llik = 
  #     (1-alpha) * (
  #       mvtnorm::dmvnorm(y, rep(0, n), cov_mat, log = TRUE) +
  #       sum(dbinom(z, 1, p_incl, log = TRUE)) -
  #       sum(dbinom(z, 1, sigmoid_phi, log = TRUE))
  #     )
  #   
  #   return(exp(llik))
  # }
  # 
  while (!converged & iter <= max_iter) {
  
    sigmoid_phi = plogis(phi)
    
    # browser()
    # log of MC estimate of alpha divergence
    sum_log_alpha_div = 0
    # sum_alpha_div = 0
    # min_log_ad = Inf
    for (k in 1:n_div) {
      z = rbinom(p, size = 1, sigmoid_phi)
      log_ad = log_alpha_div(z, phi, X, y, v_noise, v_slab, p_incl)
      
      # small value truncation
      # if (log_ad < min_logval)
      #   log_ad = min_logval
      
      # remove minimum from sum
      # if (log_ad < min_log_ad) {
      #   min_log_ad = log_ad
      # }
      
      sum_log_alpha_div = sum_log_alpha_div + log_ad
      # sum_alpha_div = sum_alpha_div + exp(log_ad)
    }
    # exp_alpha_div = sum_alpha_div / n_div
    # alpha_div_history[iter] = log(exp_alpha_div)
    alpha_div_history[iter] = sum_log_alpha_div / (1-alpha) 
    
    # browser()
    # u2g estimate of gradient
    u2g_grad = 0    
    for (k in 1:n_u2g) {
      u = runif(p)
      umask = u > (1 - sigmoid_phi)
      dmask = u < sigmoid_phi
      
      log_fu = log_alpha_div(umask, phi, X, y, v_noise, v_slab, p_incl)
      if (log_fu < min_logval) log_fu = min_logval
      log_fd = log_alpha_div(dmask, phi, X, y, v_noise, v_slab, p_incl)
      if (log_fd < min_logval) log_fd = min_logval
      fu = exp(log_fu)
      fd = exp(log_fd)

      # arm grad
      # u2g_grad = u2g_grad + (fu - fd) * (u - .5) * abs(umask - dmask)  
      
      # u2g grad
      u2g_grad = u2g_grad + .5*(fu - fd) * (umask - dmask) * max(phi, 1 - phi)
    }
    u2g_grad = - u2g_grad / n_u2g
    
    # browser()
    
    phi_old = phi
    phi = phi - lr * u2g_grad / ( (1-alpha) * exp_alpha_div )
    dphi = phi - phi_old
    
    iter = iter + 1
    
    if ( max(abs(dphi)) < delta )
      converged = TRUE
  }
  
  
  result = list(
    alpha_div = alpha_div_history,
    logodds = phi,
    p = plogis(phi),
    iters = iter
  )
  
  return(result)
}

#' @param z sampled variable inclusion mask (p x 1)
#' @param phi variational variable inclusion parameter (p x 1) 
elbo = function(z, phi, X, y, v_noise, v_slab, p_incl) {
  p = ncol(X)
  n = nrow(X)
  
  X_select = X[, z]
  cov_mat = v_slab * tcrossprod(X_select) + v_noise * diag(n)
  
  mvtnorm::dmvnorm(y, rep(0, n), cov_mat, log = TRUE) +
    sum(dbinom(z, 1, p_incl, log = TRUE)) -
    sum(dbinom(z, 1, plogis(phi), log = TRUE))
}

l0_ss = function(X, y, v_noise, v_slab, p_incl, 
                 n_u2g = 50, n_elbo = NULL, 
                 lr = .1, delta = 1e-4, max_iter = 200) {
  p = ncol(X)
  n = nrow(X)
  phi = rep(0, p)
  
  converged = FALSE
  iter = 0
  elbo_history = rep(0, max_iter)
  
  while (!converged & iter < max_iter) {
    
    sigmoid_phi = plogis(phi)
    
    # f = function(z) {
    #   cov_mat = X %*% ( (v_slab*diag(p)) * tcrossprod(z) ) %*% t(X) + v_noise*diag(n)
    #   
    #   mvtnorm::dmvnorm(y, rep(0, n), cov_mat, log = TRUE) +
    #     sum(dbinom(z, 1, p_incl, log = TRUE)) -
    #     sum(dbinom(z, 1, sigmoid_phi, log = TRUE))
    # }
    
    u2g_grad = 0
    # 
    # u = matrix(runif(n_u2g * p), nrow = p, ncol = n_u2g)
    # 
    # umask = u > plogis(-phi)
    # dmask = u < plogis(phi)
    
    for (k in 1:n_u2g) {
      u = runif(p)
      umask = u > (1 - sigmoid_phi)
      dmask = u < sigmoid_phi
      
      fu = elbo(umask, phi, X, y, v_noise, v_slab, p_incl)
      fd = elbo(dmask, phi, X, y, v_noise, v_slab, p_incl)
      
      # u2g grad
      # u2g_grad = u2g_grad + .5*(fu - fd) * plogis(abs(phi)) * (umask - dmask)
      u2g_grad = u2g_grad + .5*(fu - fd) * (umask - dmask) * max(phi, 1 - phi)
    }
    u2g_grad = - u2g_grad / n_u2g
    
    phi_old = phi
    phi = phi - lr * u2g_grad
    dphi = phi - phi_old
    
    iter = iter + 1
    if (!is.null(n_elbo)) {
      
      for (k in 1:n_elbo) {
        z = rbinom(n_elbo, size = 1, prob = plogis(phi))
        elbo_history[iter] = 
          elbo_history[iter] + elbo(z, phi, X, y, v_noise, v_slab, p_incl)
      }
      
      elbo_history[iter] = elbo_history[iter] / n_elbo
    }
    
    if ( max(abs(dphi)) < delta )
      converged = TRUE
  }
  
  result = list(
    logodds = phi,
    p = plogis(phi),
    iters = iter,
    elbo = elbo_history
  )
  
  return(result)
}