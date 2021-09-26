#' expectation propagation for Bayesian variable selection
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
epvbs = function(X, y, sigma, sa, logodds, v_inf = 100, max_iter = 200, 
                      delta = 1e-4, k = .99, damping = .5, woodbury = FALSE, opt = TRUE,
                      method = "ep2", verbose = FALSE) {
  
  p = ncol(X)
  n = nrow(X)
  ns = length(logodds)
  
  mliks = rep(NA, ns)
  post_incls = matrix(NA, p, ns)
  mu_mat = matrix(NA, p, ns)
  v_mat = matrix(NA, p, ns)
  
  sigma_vec = rep(NA, ns)
  sa_vec = rep(NA, ns)
  iters = rep(NA, ns)
  
  for (i in 1:ns) {
    v_noise = sigma[i]
    v_slab = sa[i]
    p_incl = plogis(logodds[i])
    
    if (method == "ss_ep2") {
      fit = ep_ss2(X, y, v_noise, v_slab, p_incl,
                   v_inf = v_inf,
                   max_iter = max_iter,
                   delta = delta,
                   damping = damping,
                   k = k,
                   woodbury = woodbury,
                   opt = opt)
    }
    else if (method == "ss_ep1") {
      fit = ep_ss1(X, y, v_noise, v_slab, p_incl,
                   v_inf = v_inf,
                   max_iter = max_iter,
                   delta = delta,
                   damping = damping,
                   k = k,
                   woodbury = woodbury,
                   opt = opt)
    }
    else if (method == "gss") {
      fit = GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p_incl, ncol(X)),
                              v1 = v_slab, verbose = FALSE, opt = opt,
                              damping = damping, k = k)
    }
    else {
      stop("Invalid inference method")
    }
    
    if (method %in% c("ss_ep2", "ss_ep1")) {
      mliks[i] = fit$llik
      post_incls[, i] = fit$p
      mu_mat[, i] = fit$m
      v_mat[, i] = fit$v  
    }
    else if (method == "gss") {
      mliks[i] = fit$evidence
      post_incls[, i] = plogis(fit$posteriorApproximation$p)
      mu_mat[, i] = fit$meanMarginals
      v_mat[, i] = fit$varMarginals
    }
    
    sigma_vec[i] = fit$v_noise
    sa_vec[i] = fit$v_slab
    iters[i] = fit$iters
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
    mliks = mliks,
    iters = iters
  )
  
  return(result)
}

#' group spike and slab prior
#' @export
ep_grid_gss = function(X, y, sigma, sa, logodds, verbose=FALSE, opt = TRUE,
                       damping = .5, k = .99, 
                       opt_lower = c(0, 0), opt_upper = c(Inf, Inf)) {
  
  p = ncol(X)
  n = nrow(X)
  ns = length(logodds)
  
  mliks = rep(NA, ns)
  logodds_incls = matrix(NA, p, ns)
  mu_mat = matrix(NA, p, ns)
  v_mat = matrix(NA, p, ns)
  
  sigma_vec = rep(NA, ns)
  sa_vec = rep(NA, ns)
  
  for (i in 1:ns) {
    
    if (verbose) {
      print(paste0("Hyperparameter setting ", i))
    }
    
    v_noise = sigma[i]
    v_slab = sa[i]
    p_incl = plogis(logodds[i])
    
    fit = GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p_incl, ncol(X)),
                            v1 = v_slab, verbose = FALSE, opt = opt,
                            damping = damping, k = k)
    
    mliks[i] = fit$evidence
    logodds_incls[, i] = fit$posteriorApproximation$p
    mu_mat[, i] = fit$meanMarginals
    v_mat[, i] = fit$varMarginals
    
    sigma_vec[i] = fit$v_noise
    sa_vec[i] = fit$v_slab
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

#' @export
vb_grid_ss = function(X, y, sigma, sa, logodds, verbose=FALSE, opt=TRUE) {
  
  p = ncol(X)
  n = nrow(X)
  ns = length(logodds)
  
  elbos = rep(NA, ns)
  logodds_incls = matrix(NA, p, ns)
  mu_mat = matrix(NA, p, ns)
  v_mat = matrix(NA, p, ns)
  n_iters = rep(NA, ns)
  
  sigma_vec = rep(NA, ns)
  sa_vec = rep(NA, ns)
  
  for (i in 1:ns) {
    
    if (verbose) {
      print(paste0("Hyperparameter setting ", i))
    }
    
    v_noise = sigma[i]
    v_slab = sa[i]
    p_incl = plogis(logodds[i])
    
    fit = vb_ss(X, y, v_noise, v_slab, p_incl, opt=opt)
    
    elbos[i] = fit$elbo
    n_iters[i] = fit$iters
    logodds_incls[, i] = qlogis(fit$alpha)
    mu_mat[, i] = fit$mu
    v_mat[, i] = fit$v
    
    sigma_vec[i] = fit$v_noise
    sa_vec[i] = fit$v_slab
  }
  
  post_incls = plogis(logodds_incls)
  
  weights = normalizelogweights(elbos)
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
    mliks = elbos,
    iters = n_iters
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