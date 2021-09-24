#' Batch variational bayes for linear regression spike and slab prior
#' 
#' Notice the updates are the same as in Carbonetto and Stephens 2012, only 
#' difference is batchwise order of updates. 
#' 
#' @param X n x p matrix of covariates
#' @param y n x 1 vector of response
#' @param v_noise scalar observation variance hyperparameter
#' @param v_slab scalar slab variance hyperparameter
#' @param p_incl scalar or vector prior covariate inclusion probability hyperparameter
#' @param delta convergence threshold
#' @param eps perturbation to avoid NaN for 0 or 1 posterior inclusion probability
#' @export
vb_ss = function(X, y, v_noise, v_slab, p_incl, max_iter = 200, delta = 1e-4,
                 prob_thresh = 1e-7, eps = 1e-6, opt=TRUE) {
  
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
    
    # updates (E step)
    v = v_noise / ((1 / v_slab) + diag_XtX)
    
    # browser()
    mu = (v / v_noise) * (Xty - XtX %*% (mu * alpha) + diag_XtX * (mu * alpha))
    
    logit_alpha = qlogis(p_incl) + .5 * (mu^2 / v) + .5*log(v / (v_noise * v_slab))
    # too_small = logit_alpha < qlogis(prob_thresh)
    # too_large = logit_alpha > qlogis(1 - prob_thresh)
    alpha = plogis(logit_alpha)
    # alpha[too_small] = 
    
    # hyperparameter optimization (M step)
    mu_alpha = mu * alpha
    sse = sum( (y - X %*% (mu_alpha))^2 )
    alpha_vmu2 = alpha * (v + mu^2)
    var_beta = diag_XtX %*% (alpha_vmu2 - mu_alpha^2)
    if (opt) {
      v_noise = c(( sse + var_beta + sum(alpha_vmu2) / v_slab ) / (n + sum(alpha) ))
      v_slab = c( sum(alpha_vmu2) / (v_noise * sum(alpha)) )
    }
    
    # elbo
    e1 = -.5*log(2*pi*v_noise)
    e2 = - sse / (2*v_noise)
    e3 = - var_beta / (2*v_noise)
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
    v = v,
    v_noise = v_noise,
    v_slab = v_slab
  )
  
  return(result)
}