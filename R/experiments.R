discrete_independent = function(n_sim, output_path,
                                n_obs, covariates, covariate_dim, tau,
                                lambda, nugget = 10,
                                p_incl_grid, v_noise_init = 1, v_slab_init = 3, 
                                seed = 1) {
  n_cov = length(covariates)
  p = length(lambda) - 1
  
  Var = solve(lambda %*% t(lambda) + diag(rep(nugget, p+1)))
  
  covariates_full = matrix(rep(covariates, each = n_obs / n_cov),
                           nrow = n_obs, ncol = covariate_dim, byrow = FALSE)
  covariates_levels = matrix(covariates, nrow = n_cov, ncol = covariate_dim, byrow = FALSE)
  weight_idx = seq(1, n_obs * (n_cov-1) / n_cov + 1, n_obs / n_cov)
  weight_mat = epwpl::weight_matrix(covariates_full, tau)[weight_idx, ]
  
  # true graph
  true_graph = matrix(0, p+1, p+1)
  for(i in 1:(p+1)){
    for(j in 1:(p+1)){
      true_graph[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
    }}
  diag(true_graph) = 0
  true_individual_graphs = replicate(2, true_graph, simplify = FALSE)
  
  # hyperparameters
  v_noise = v_noise_init
  v_slab = v_slab_init
  
  set.seed(1)
  progressr::with_progress({
    prog = progressr::progressor(along = 1:n_sim)
    
    metrics = foreach::foreach(sim_idx = 1:n_sim, .combine = rbind) %dorng% {
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate the data for this iteration
      data_mat = MASS::mvrnorm(n_obs, rep(0, p+1), Var)
      
      varbvs_vopt_result =
        epwpl::wpl_varbvs(data_mat,
                          covariates = covariates_levels, tau = tau, weight_mat = weight_mat,
                          v_noise_grid = v_noise,
                          v_slab_grid = v_slab,
                          p_incl_grid = p_incl_grid,
                          opt = TRUE)
      
      ep_vopt_result =
        epwpl::wpl_ep(data_mat,
                      covariates = covariates_levels, tau = tau, weight_mat = weight_mat,
                      v_noise_grid = v_noise,
                      v_slab_grid = v_slab,
                      p_incl_grid = p_incl_grid,
                      damping = .9, k = .99,
                      opt = TRUE,
                      opt_method = "Nelder-Mead")
    }
    
    metrics = rbind(
      epwpl::score_graphs(ep_vopt_result, true_individual_graphs),
      epwpl::score_graphs(varbvs_vopt_result, true_individual_graphs),
      epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .25, true_individual_graphs),
      epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .5, true_individual_graphs),
      epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .75, true_individual_graphs)
    )
    metrics[, sim_idx := sim_idx]
  })
}