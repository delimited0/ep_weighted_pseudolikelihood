#' @export
wpl_ep = function(data_mat, weight_mat, 
                  v_noise_grid, v_slab_grid, p_incl_grid, 
                  damping = .5, k = .99,
                  opt = TRUE, verbose = FALSE) {
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
      
      fit = epwpl::ep_grid_ss(X_weighted, y_weighted, 
                               v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                               eps = damping, k = k,
                               opt = opt)
      
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
  
  n = nrow(data_mat)
  p = ncol(data_mat)-1
  
  mylist = rep(list(matrix(0, n, p)), p+1)
  
  dim_elbos = rep(NA, p+1)
  
  for (resp_index in 1:(p+1)) {
    
    y = data_mat[, resp_index] #Set variable number `resp_index` as the response
    X_mat = data_mat[, -resp_index]
    X_vec = matrix(0, n*p, 1)
    X = matrix(rep(0, n^2*p), nrow = n, ncol = n*p)
    
    for(i in 1:n) {
      for(j in 1:p) {
        k = p*(i-1) + 1 
        X[i, k+j-1] = X_mat[i, j]
        X_vec[k+j-1] = X[i, k+j-1]
      }
    }
    
    sigmasq = sigma0^2
    
    XtX = t(X) %*% X
    
    DXtX = diag(XtX)
    DXtX_rep = rep(DXtX, p) 
    DXtX_mat = matrix(DXtX_rep, n*p, p, byrow=FALSE)
    
    alpha = rep(.2, n*p)
    sigmabeta_sq = v_slab
    mu = rep(0, p) # Variational parameter
    true_pi = p0
    
    y_long_vec = as.vector(t(y %*% matrix(1, 1, p)))
    Xty=t(X) %*% as.vector(y)
    mu_mat = matrix(NA, n, p)
    
    D_long = matrix(0, n*p, n)
    for( i in 1:n){
      D_long[, i] = matrix(t(weight_mat[, i] %*% matrix(1, 1, p)), n*p, 1)
    }
    
    ind_vec = seq(0,(n-1)*p,by=p)
    Ind_mat = matrix(0,n,p)
    for(j in 1:p){
      Ind_mat[, j] = ind_vec + j
    }
    Big_ind = matrix(0, n*p, p)
    Big_ind_1 = matrix(0, n*p, p)
    for(j in 1:p){
      Big_ind[Ind_mat[, j], j] = 1
      Big_ind_1[Ind_mat[, j], j] = 0
    }
    
    DXtX_Big_ind = DXtX_mat * Big_ind
    
    candL = seq(0.1, 0.9, .2)# Different values of hyperparameter true_pi
    #candL=0.5
    like = rep(0, length(candL))
    elb = like
    
    est_pi = rep(0, n)
    est_q = est_pi
    
    
  }
  
    
}

