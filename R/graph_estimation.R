#' @export
wpl_ep = function(data_mat, covariates, tau,
                  v_noise_grid, v_slab_grid, p_incl_grid, 
                  damping = .5, k = .99,
                  woodbury = FALSE, 
                  opt = TRUE, opt_method = "Nelder-Mead",
                  lb = c(0, 0), ub = c(Inf, Inf),
                  verbose = FALSE, method = "ss_ep2",
                  n_threads = 1) {
  # registerDoFuture()
  # plan(multisession, workers = n_threads)
  # RhpcBLASctl::blas_set_num_threads(blas_threads)
  # RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  n = nrow(covariates)
  
  weight_mat = weight_matrix(covariates, tau)
  sqrt_weight = sqrt(weight_mat)
  
  individuals = vector("list", n)

  model_llik = 0
  
  # progressr::with_progress({
  #   prog = progressr::progressor(along = 1:n)
  
  for (i in 1:n) {
    # prog(sprintf("Individual =%g", i))
    
    graph = matrix(0, p, p)
    fits = vector("list", p)
    llik = 0
    
    for (resp_idx in 1:p) {
      
      print(sprintf("Individual: %d, dimension: %d", i, resp_idx))
      
      y = data_mat[, resp_idx]
      X = data_mat[, -resp_idx]
      y_weighted = y * sqrt_weight[i, ]
      X_weighted = X * sqrt_weight[i, ]
      
      fit = epwpl::epvbs(X_weighted, y_weighted, 
                         v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                         damping = damping, k = k,
                         woodbury = woodbury, opt = opt, opt_method = opt_method,
                         method = method)
      
      graph[resp_idx, -resp_idx] = t(fit$pip)
      fits[[resp_idx]] = fit
      llik = llik + fit$llik
    }
    
    individuals[[i]]$graph = graph
    individuals[[i]]$fits = fits
    individuals[[i]]$llik = llik
    individuals[[i]]$covariates = covariates[i, ]
    
    model_llik = model_llik + llik
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
         method_name = method,
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
wpl_varbvs = function(data_mat, covariates, tau,
                      v_noise_grid, v_slab_grid, p_incl_grid,
                      opt = TRUE, verbose = FALSE) {
  
  p = ncol(data_mat)
  n = nrow(covariates)
  
  weight_mat = weight_matrix(covariates, tau)
  sqrt_weight = sqrt(weight_mat)
  
  individuals = vector("list", n)
  
  model_llik = 0
  
  for (i in 1:n) {
    
    graph = matrix(0, p, p)
    fits = vector("list", p)
    llik = 0
    
    for (resp_idx in 1:p) {
      
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
      llik = llik + fit$llik
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

