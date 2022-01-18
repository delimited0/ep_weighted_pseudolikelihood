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
    
    param_names = list("target_dim" = 1:p, "predictor_dim" = 1:(p-1), "hyp_idx" = 1:ns)
    alpha = array(NA, dim = c(p, p-1, ns), dimnames = param_names)
    mu = array(NA, dim = c(p, p-1, ns), dimnames = param_names)
    v = array(NA, dim = c(p, p-1, ns), dimnames = param_names)
    # v_site2 = array(NA, dim = c(p, p-1, ns), dimnames = param_names)
    Cov = array(NA, dim = c(p, p-1, p-1, ns), 
                dimnames = list("target_dim" = 1:p, "predictor_dim1" = 1:(p-1),
                                "predictor_dim2" = 1:(p-1), "hyp_idx" = 1:ns))
    hyp = array(NA, dim = c(p, 3, ns), 
                dimnames = list("target_dim" = 1:p,
                                "hyp_param" = c("v_noise", "v_slab", "logodds"),
                                "hyp_idx" = 1:ns))
    
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
      
      alpha[resp_idx, , ] = fit$alpha
      mu[resp_idx, , ] = fit$mu
      v[resp_idx, , ] = fit$v
      # v_site2[resp_idx, , ] = fit$v_site2
      Cov[resp_idx, , , ] = fit$Cov
      hyp[resp_idx, , ] = rbind(fit$sigma, fit$sa, fit$logodds)
    }
    
    individuals[[i]]$graph = graph
    individuals[[i]]$fits = fits
    individuals[[i]]$llik = llik
    individuals[[i]]$alpha = alpha
    individuals[[i]]$mu = mu
    individuals[[i]]$v = v
    # individuals[[i]]$v_site2 = v_site2
    individuals[[i]]$Cov = Cov
    individuals[[i]]$hyp = hyp
    
    individuals[[i]]$covariates = covariates[i, ]
  }
  
  result = list(individuals = individuals,
                tau = tau,
                n = n, p = p,
                opt = opt,
                opt_method = opt_method,
                damping = damping,
                k = k,
                method_name = method)
  
  return(structure(result, class = c("ep_wpl", "wpl")))
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

sample_posterior = function(graph_est) UseMethod("sample_posterior")

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
        mu_given_z = fit$mu[, hyp_idx]
        mu_given_z[z == 0] = 0
        s_given_z = fit$s[, hyp_idx]
        s_given_z[z == 0] = 0
        
        w = rnorm(p-1, mu_given_z, s_given_z)
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

sample_ep_posterior = function(n_samples, ep_graph_est) {
  
  n = length(ep_graph_est$individuals)
  p = ep_graph_est$p
  
  dnames = list("sample" = 1:n_samples,
                "individual" = 1:n,
                "target_dim" = 1:p,
                "predictor_dim" = 1:(p-1))
  
  hyp_samples = 
    array(NA, dim = c(n_samples, n, p, 3), 
          dimnames = c(dnames[-4], list("hyp_param" = c("v_noise", "v_slab", "p_incl")))
    )
  gamma_samples = array(NA, dim = c(n_samples, n, p, p-1), dimnames = dnames)
  beta_samples = array(NA, dim = c(n_samples, n, p, p-1), dimnames = dnames)
  
  for (i in 1:n) {
    
    individual = ep_graph_est$individuals[[i]]
    
    for (j in 1:p) {
      
      # fit = individual$fits[[j]]
      hyp_prob = epwpl::normalizelogweights(individual$llik[j, ])
      ns = length(hyp_prob)
      
      for (k in 1:n_samples) {
        
        hyp_idx = sample(x = 1:ns, size = 1, prob = hyp_prob)
        z = rbinom(p-1, 1, individual$alpha[j, , hyp_idx])
        w = MASS::mvrnorm(n = 1, individual$mu[j, , hyp_idx], individual$Cov[j, , , hyp_idx])  
        
        hyp_samples[k, i, , ] = individual$hyp[, , hyp_idx]
        gamma_samples[k, i, j, ] = z
        beta_samples[k, i, j, ] = w
      }
    }
  }
  
  return(list(
    gamma = gamma_samples,
    beta = beta_samples,
    hyp = hyp_samples
  ))
}

#' @export 
log_joint_prob = function(X, weight_mat, beta_samples, gamma_samples, 
                          v_noise, v_slab) {
  
  n = nrow(weight_mat)
  
  sqrt_weight = sqrt(weight_mat)
  p = ncol(X)
  
  for (l in 1:n) {
    
    for (tgt_idx in 1:p) {
      
      target = X[, tgt_idx] * sqrt_weight[l, ]
      preds = X[, -tgt_idx] * sqrt_weight[l, ]
      beta_preds = beta[-tgt_idx]
      
      dnorm(target, preds %*% beta_preds, )
            
    }
  }
  
  
}
