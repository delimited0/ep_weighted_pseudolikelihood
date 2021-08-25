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
                    delta = 1e-4, k = .99, woodbury = FALSE) {
  p = ncol(X)
  n = nrow(X)
  ns = length(sigma)
  
  mliks = rep(NA, p)
  post_incls = matrix(NA, p, ns)
  alpha_mat = matrix(NA, p, ns)
  
  for (i in 1:ns) {
    sigma0 = sigma[i]
    s_slab = sa[i]
    p0 = plogis(logodds[i])
    
    fit = ep_wlr(X, y, sigma0, p0, s_slab, 
                 v_inf = v_inf, 
                 max_iter = max_iter,
                 delta = delta,
                 k = k,  
                 woodbury = woodbury,
                 opt = FALSE)
    
    mliks[i] = fit$llik
    post_incls[, i] = fit$p
    alpha_mat[, i] = fit$p
  }
  
  weights = mliks / sum(mliks)
  
  pip = post_incls * weights
  
  result = list(
    pip = pip,
    w = weights,
    
  )
  
}