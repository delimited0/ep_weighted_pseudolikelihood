#' @export
vb_ss = function(X, y, v_noise, v_slab, p_incl, max_iter = 200, delta = 1e-4,
                 prob_thresh = 1e-7, eps = 1e-6) {
  
  p = ncol(X)
  n = length(y)
  
  mu = rep(0, p)
  v = rep(0, p)
  logit_alpha = rep(0, p)
  alpha = rep(p_incl, p)
  
  converged = FALSE
  iter = 0
  elbos = rep(NA, max_iter)
  
  Xty = t(X) %*% y
  XtX = t(X) %*% X
  diag_XtX = diag(XtX)
  
  while (!converged & iter <= max_iter) {
    
    iter = iter + 1
    alpha_old = alpha
    
    # updates
    v = v_noise / ((1 / v_slab) + diag_XtX)
    
    # browser()
    mu = (v / v_noise) * (Xty - XtX %*% (mu * alpha) + diag_XtX * (mu * alpha))
    
    logit_alpha = qlogis(p_incl) + .5 * (mu^2 / v) + .5*log(v / (v_noise * v_slab))
    # too_small = logit_alpha < qlogis(prob_thresh)
    # too_large = logit_alpha > qlogis(1 - prob_thresh)
    alpha = plogis(logit_alpha)
    # alpha[too_small] = 
    
    mu_alpha = mu * alpha
    e1 = -.5*log(2*pi*v_noise)
    e2 = - sum( (y - X %*% (mu_alpha))^2 ) / (2*v_noise)
    e3 = - diag_XtX %*% (alpha * (v + mu^2) - mu_alpha^2) / (2*v_noise)
    e4 = - sum(alpha * (log(alpha + eps) - log(p_incl)))
    e5 = - sum((1 - alpha) * (log(1 - alpha + eps) - log(1 - p_incl)))
    e6 = .5*sum(alpha * (1 + log(v) - log(v_slab * v_noise) - 
                           (v + mu^2) / (v_slab * v_noise)))
    
    elbo = e1 + e2 + e3 + e4 + e5 + e6
    elbos[iter] = elbo
    
    # browser()
    change_alpha = alpha - alpha_old
    if (all(abs(change_alpha) < delta))
      converged = TRUE
  }
  
  result = list(
    elbo = elbo,
    elbo_history = elbos[!is.na(elbos)],
    iters = iter,
    mu = mu,
    alpha = alpha,
    v = v
  )
  
  return(result)
}