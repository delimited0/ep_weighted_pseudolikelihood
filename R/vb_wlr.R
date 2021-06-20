
#" batch VB as in Huang et al 2016
#' @export
vb_wlr = function(X, y, sigma_noise, p0, v_slab, 
                  mu_init = rep(0, p), 
                  max_iter = 200, delta = 1e-4, thresh = 1e-7) {
  
  p = ncol(X)
  
  mu = mu_init
  alpha = rep(p0, p)
  logit_alpha = rep(qlogis(p0), p)
  
  logit_p0 = qlogis(p0) 
  lthresh = qlogis(thresh)
  uthresh = qlogis(1 - thresh)
  
  # precompute
  # variational variance 
  ssq = sigma_noise^2 / ((1 / v_slab) + colSums(X^2))
  
  XtX = t(X) %*% X
  Xty = t(X) %*% y
  
  converged = FALSE
  iter = 1
  
  while (!converged & iter < max_iter) {
    
    mu_old = mu
    logit_alpha_old = logit_alpha
    alpha_old = alpha
    
    mu = (ssq / sigma_noise^2) * (Xty - XtX %*% (mu * alpha))
    
    logit_alpha = 
      logit_p0 + 
      .5 * (mu^2 / ssq) +
      log( sqrt(ssq) / (sigma_noise * sqrt(v_slab)) )
    
    # threshhold alpha for stability
    # idx_large = logit_alpha > uthresh
    # idx_small = logit_alpha < lthresh
    # logit_alpha[idx_large] = uthresh
    # logit_alpha[idx_small] = lthresh
    # 
    # alpha[logit_alpha > 9] = 1
    # alpha[logit_alpha <= 9] = plogis(logit_alpha[logit_alpha <= 9])
    
    # elbo_opt = 
    #   dfoptim::nmkb(c(sigma0, v_slab), hyper_obj,
    #                 lower = c(0, 0), upper = c(Inf, Inf),
    #                 X = X, y = y, ssq = ssq, mu = mu, alpha = alpha, p0 = p0)
    # 
    # sigma_noise = elbo_opt$par[1]
    # v_slab = elbo_opt$par[2]
    
    # check convergence
    if (all(abs(mu - mu_old) < delta) &
        all(abs(logit_alpha - logit_alpha_old) < delta))
    # if (sqrt(sum((alpha - alpha_old)^2)) < delta)
      converged = TRUE
    
    iter = iter + 1
  }
  
  result <- list(
    mu = mu,
    logit_alpha = logit_alpha,
    var = ssq,
    iters = iter,
    sigma0 = sigma0,
    v_slab = v_slab
  )
}

# elbo = function(X, y, ssq, mu, alpha, sigma_noise, p0, v_slab) {
#   browser()
#   t1 = -sum( (y - X %*% (mu * alpha))^2 ) / (2*sigma_noise^2)
#   t2 = -sum( X^2 %*% (alpha * (mu^2 + ssq) - alpha^2 * mu^2) ) / (2*sigma_noise^2)
#   t3 = sum( alpha * ((1 + log(ssq))) ) / 2
#   t4 = -sum( alpha * log((alpha + 1e-5) / p0) + 
#                (1-alpha)*log((1-alpha + 1e-5) / (1 - p0)) )
#   t5 = -sum( alpha * ((mu^2 + ssq) / (2*sigma_noise^2*v_slab) + 
#                         log(sigma_noise^2*v_slab) / 2) )
#   t6 = sum( .5*log(1 / (2 * pi * sigma_noise^2)) )
# 
#   return(t1 + t2 + t3 + t4 + t5 + t6)
# }

elbo_obj = function(params, X, y, ssq, mu, alpha, p0) {
  sigma_noise = params[1]
  v_slab = params[2]
  
  t1 = -sum( (y - X %*% (mu * alpha))^2 ) / (2*sigma_noise^2)
  t2 = -sum( X^2 %*% (alpha * (mu^2 + ssq) - alpha^2 * mu^2) ) / (2*sigma_noise^2)
  t3 = sum( alpha * ((1 + log(ssq))) ) / 2
  t4 = -sum( alpha * log((alpha + 1e-5) / p0) + 
               (1-alpha)*log((1-alpha + 1e-5) / (1 - p0)) )
  t5 = -sum( alpha * ((mu^2 + ssq) / (2*sigma_noise^2*v_slab) + 
                        log(sigma_noise^2*v_slab) / 2) )
  t6 = sum( .5*log(1 / (2 * pi * sigma_noise^2)) )
  
  return(t1 + t2 + t3 + t4 + t5 + t6)
}