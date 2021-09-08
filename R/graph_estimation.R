#' @export
wpl_ep = function(data_mat, weight_mat, 
                  v_noise_grid, v_slab_grid, p_incl_grid, 
                  opt = TRUE) {
  # registerDoFuture()
  # plan(multisession, workers = n_threads)
  # RhpcBLASctl::blas_set_num_threads(blas_threads)
  # RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  n = nrow(weight_mat)
  
  sqrt_weight = sqrt(weight_mat)
  
  graphs = replicate(n, matrix(0, p, p), simplify = FALSE)
  
  # llik = 0
  
  progressr::with_progress({
    prog = progressr::progressor(along = 1:n)
    
    for (i in 1:n) {
      prog(sprintf("Individual =%g", i))
      
      for (resp_idx in 1:p) {
        
        y = data_mat[, resp_idx]
        X = data_mat[, -resp_idx]
        y_weighted = y * sqrt_weight[i, ]
        X_weighted = X * sqrt_weight[i, ]
        
        fit = epwpl::ep_grid_gss(X_weighted, y_weighted, 
                                 v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                                 verbose=FALSE)
        
        graphs[[i]][resp_idx, -resp_idx] = t(plogis(fit$pip))
        
        
        # llik = llik + fit$llik
      }
    }
  })
  
  return(list(graphs = graphs
              # llik = llik
              ))
}