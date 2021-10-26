#' @export
wpl_ep = function(data_mat, covariates, tau, weight_mat = NULL,
                  v_noise_grid, v_slab_grid, p_incl_grid, 
                  damping = .5, k = .99,
                  woodbury = FALSE, 
                  opt = TRUE, opt_method = "Nelder-Mead",
                  lb = c(0, 0), ub = c(Inf, Inf),
                  verbose = FALSE, method = "ss_ep2",
                  n_threads = 1, blas_threads = 1) {
  # registerDoFuture()
  # plan(multisession, workers = n_threads)
  RhpcBLASctl::blas_set_num_threads(blas_threads)
  RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  ns = length(p_incl_grid)
  
  if (is.null(weight_mat)) 
    weight_mat = weight_matrix(covariates, tau)
  sqrt_weight = sqrt(weight_mat)
  n = nrow(sqrt_weight)
  
  individuals = vector("list", n)

  # progressr::with_progress({
  #   prog = progressr::progressor(along = 1:n)
  
  for (i in 1:n) {
    # prog(sprintf("Individual =%g", i))
    
    graph = matrix(0, p, p)
    fits = vector("list", p)
    llik = matrix(NA, p, ns)
    alpha = matrix(NA, p, ns)
    mu = matrix(NA, p, ns)
    
    
    for (resp_idx in 1:p) {
      
      if (verbose)
        print(sprintf("Individual: %d, dimension: %d", i, resp_idx))
      
      y = data_mat[, resp_idx]
      X = data_mat[, -resp_idx]
      y_weighted = y * sqrt_weight[i, ]
      X_weighted = X * sqrt_weight[i, ]
      
      fit = epwpl::epbvs(X_weighted, y_weighted, 
                         v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                         damping = damping, k = k,
                         woodbury = woodbury, opt = opt, opt_method = opt_method,
                         method = method)
      
      graph[resp_idx, -resp_idx] = t(fit$pip)
      fits[[resp_idx]] = fit
      llik[resp_idx, ] = fit$mliks
    }
    
    individuals[[i]]$graph = graph
    individuals[[i]]$fits = fits
    individuals[[i]]$llik = llik
    individuals[[i]]$covariates = covariates[i, ]
  }
  # })
  
  return(
    list(individuals = individuals,
         tau = tau,
         n = n, p = p,
         opt = opt,
         opt_method = opt_method,
         damping = damping,
         k = k,
         method_name = method
    )
  )
}


#' @export
wpl_varbvs = function(data_mat, covariates, tau, weight_mat = NULL,
                      v_noise_grid, v_slab_grid, p_incl_grid,
                      opt = TRUE, verbose = FALSE, blas_threads = 1) {
  RhpcBLASctl::blas_set_num_threads(blas_threads)
  RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  ns = length(p_incl_grid)
  
  if (is.null(weight_mat)) 
    weight_mat = weight_matrix(covariates, tau)
  sqrt_weight = sqrt(weight_mat)
  n = nrow(sqrt_weight)
  
  individuals = vector("list", n)
  
  model_llik = 0
  
  for (i in 1:n) {
    
    graph = matrix(0, p, p)
    fits = vector("list", p)
    llik = matrix(NA, p, ns)
    
    for (resp_idx in 1:p) {
      
      if (verbose)
        print(sprintf("Individual: %d, dimension: %d", i, resp_idx))
      
      y = data_mat[, resp_idx]
      X = data_mat[, -resp_idx]
      y_weighted = y * sqrt_weight[i, ]
      X_weighted = X * sqrt_weight[i, ]
    
      fit = varbvs::varbvs(X = X, Z = NULL, y = y, 
                           weights = as.vector(sqrt_weight[i, ]),
                           sigma = v_noise_grid,
                           sa = v_slab_grid,
                           logodds = qlogis(p_incl_grid),
                           update.sigma = opt,
                           update.sa = opt,
                           verbose = verbose)
      
      graph[resp_idx, -resp_idx] = t(fit$pip)
      fits[[resp_idx]] = fit
      llik[resp_idx, ] = fit$logw
    }
    
    individuals[[i]]$graph = graph
    individuals[[i]]$fits = fits
    individuals[[i]]$llik = llik
    individuals[[i]]$covariates = covariates[i, ]
    
    model_llik = model_llik + llik
  }
  
  return(
    list(individuals = individuals,
         tau = tau,
         n = n, p = p,
         opt = opt,
         method_name = "varbvs",
         # covariates = covariates,
         llik = model_llik
    )
  )
}

#' @export
wpl_ep_gss = function(data_mat, weight_mat, 
                      v_noise_grid, v_slab_grid, p_incl_grid, 
                      damping = .5, k = .99,
                      opt = TRUE, verbose = FALSE,
                      opt_lower = c(0, 0), opt_upper = c(Inf, Inf)) {
  # registerDoFuture()
  # plan(multisession, workers = n_threads)
  # RhpcBLASctl::blas_set_num_threads(blas_threads)
  # RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  n = nrow(weight_mat)
  
  sqrt_weight = sqrt(weight_mat)
  
  graphs = replicate(n, matrix(0, p, p), simplify = FALSE)
  results = replicate(n, vector("list", p), simplify=FALSE)
  
  # llik = 0
  
  # progressr::with_progress({
  #   prog = progressr::progressor(along = 1:n)
  
  for (i in 1:n) {
    # prog(sprintf("Individual =%g", i))
    
    for (resp_idx in 1:p) {
      
      print(sprintf("Individual: %d, dimension: %d", i, resp_idx))
      
      y = data_mat[, resp_idx]
      X = data_mat[, -resp_idx]
      y_weighted = y * sqrt_weight[i, ]
      X_weighted = X * sqrt_weight[i, ]
      
      fit = epwpl::ep_grid_gss(X_weighted, y_weighted, 
                               v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                               damping = damping, k = k,
                               verbose=verbose, opt = opt)
      
      graphs[[i]][resp_idx, -resp_idx] = t(fit$pip)
      results[[i]][[resp_idx]] = fit
      # llik = llik + fit$llik
    }
  }
  # })
  
  return(list(graphs = graphs,
              results = results
              # llik = llik
  ))
}

#' @export
wpl_vb = function(data_mat, weight_mat,
                  v_noise_grid, v_slab_grid, p_incl_grid,
                  opt = TRUE, verbose = FALSE) {
  
  p = ncol(data_mat)
  n = nrow(weight_mat)
  
  sqrt_weight = sqrt(weight_mat)
  
  graphs = replicate(n, matrix(0, p, p), simplify = FALSE)
  results = replicate(n, vector("list", p), simplify=FALSE)
  
  # llik = 0
  
  # progressr::with_progress({
  #   prog = progressr::progressor(along = 1:n)
  
  for (i in 1:n) {
    # prog(sprintf("Individual =%g", i))
    
    for (resp_idx in 1:p) {
      
      print(sprintf("Individual: %d, dimension: %d", i, resp_idx))
      
      y = data_mat[, resp_idx]
      X = data_mat[, -resp_idx]
      y_weighted = y * sqrt_weight[i, ]
      X_weighted = X * sqrt_weight[i, ]
      
      fit = epwpl::vb_grid_ss(X_weighted, y_weighted, 
                              v_noise_grid, 
                              v_slab_grid, 
                              p_incl_grid,
                              verbose=verbose, opt=opt)
      
      graphs[[i]][resp_idx, -resp_idx] = t(fit$pip)
      results[[i]][[resp_idx]] = fit
      # llik = llik + fit$llik
    }
  }
  # })
  
  return(list(graphs = graphs,
              results = results
              # llik = llik
  ))
}

#' @export
sample_vb_posterior = function(n_samples, vb_graph_est) {
  
  n = length(vb_graph_est$individuals)
  p = vb_graph_est$p
  # samples = vector("list", n)
  dnames = list("sample" = 1:n_samples,
                "individual" = 1:n,
                "target_dim" = 1:p,
                "predictor_dim" = 1:(p-1))
  gamma_samples = array(NA, dim = c(n_samples, n, p, p-1), dimnames = dnames)
  beta_samples = array(NA, dim = c(n_samples, n, p, p-1), dimnames = dnames)
  
  for (i in 1:n) {
    
    individual = vb_graph_est$individuals[[i]]
    
    for (j in 1:p) {
      
      fit = individual$fits[[j]]
      hyp_prob = epwpl::normalizelogweights(fit$logw)
      ns = length(hyp_prob)
      
      for (k in 1:n_samples) {
        
        hyp_idx = sample(x = 1:ns, size = 1, prob = hyp_prob)
        z = rbinom(p-1, 1, fit$alpha[, hyp_idx])
        w = rnorm(p-1, fit$mu[, hyp_idx], sqrt(fit$s[, hyp_idx]))
        gamma_samples[k , i, j, ] = z
        beta_samples[k, i, j, ] = w  
      }
    }
  }
  
  return(list(
    gamma = gamma_samples,
    beta = beta_samples
  ))
}

#' @export 
log_pseudo_lik = function(X, w, z) {
  
}
