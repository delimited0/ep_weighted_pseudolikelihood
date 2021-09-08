#'
#' @param X 
#' @param y
#' @param sigma
#' @param sa
#' @param logodds 
#' 
#'
#' uniform discrete grid hyperprior
#' assume prior logodds all the same for now
ep_grid_lr = function(X, y, sigma, sa, logodds, v_inf = 100, max_iter = 200, 
                      delta = 1e-4, k = .99, woodbury = FALSE, opt = TRUE) {
  
  p = ncol(X)
  n = nrow(X)
  ns = length(sigma)
  
  mliks = rep(NA, ns)
  post_incls = matrix(NA, p, ns)
  mu_mat = matrix(NA, p, ns)
  v_mat = matrix(NA, p, ns)
  
  sigma_vec = rep(NA, ns)
  sa_vec = rep(NA, ns)
  
  for (i in 1:ns) {
    v_noise = sigma[i]
    v_slab = sa[i]
    p_incl = plogis(logodds[i])
    
    fit = ep_wlr(X, y, v_noise, v_slab, p_incl, 
                 v_inf = v_inf, 
                 max_iter = max_iter,
                 delta = delta,
                 k = k,  
                 woodbury = woodbury,
                 opt = opt)
    
    mliks[i] = fit$llik
    post_incls[, i] = fit$p
    mu_mat[, i] = fit$m
    v_mat[, i] = fit$v
    
    sigma_vec[i] = fit$v_noise
    sa_vec[i] = fit$v_slab
  }
  
  weights = normalizelogweights(mliks)
  pip = post_incls %*% weights
  beta = mu_mat %*% weights
  
  result = list(
    sigma = sigma_vec,
    sa = sa_vec,
    w = weights,
    alpha = post_incls,
    mu = mu_mat,
    v = v_mat,
    pip = pip,
    beta = beta,
    weights = weights,
    mliks = mliks
  )
  
  return(result)
}

#' group spike and slab prior
#' @export
ep_grid_gss = function(X, y, sigma, sa, logodds, verbose=FALSE) {
  
  p = ncol(X)
  n = nrow(X)
  ns = length(sigma)
  
  mliks = rep(NA, ns)
  logodds_incls = matrix(NA, p, ns)
  mu_mat = matrix(NA, p, ns)
  v_mat = matrix(NA, p, ns)
  
  sigma_vec = rep(NA, ns)
  sa_vec = rep(NA, ns)
  
  for (i in 1:ns) {
    v_noise = sigma[i]
    v_slab = sa[i]
    p_incl = plogis(logodds[i])
    
    fit = GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p_incl, ncol(X)),
                            v1 = v_slab, verbose = verbose)
    
    mliks[i] = fit$evidence
    logodds_incls[, i] = fit$posteriorApproximation$p
    mu_mat[, i] = fit$meanMarginals
    v_mat[, i] = fit$varMarginals
    
    sigma_vec[i] = sigma[i]
    sa_vec[i] = sa[i]
  }
  
  post_incls = plogis(logodds_incls)
  
  weights = normalizelogweights(mliks)
  pip = post_incls %*% weights
  beta = mu_mat %*% weights
  
  result = list(
    sigma = sigma_vec,
    sa = sa_vec,
    w = weights,
    alpha = post_incls,
    mu = mu_mat,
    v = v_mat,
    pip = pip,
    beta = beta,
    weights = weights,
    mliks = mliks
  )
}

normalizelogweights = function(logw) {
  
  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  c <- max(logw)
  w <- exp(logw - c)
  
  # Normalize the probabilities.
  return(w/sum(w))
}