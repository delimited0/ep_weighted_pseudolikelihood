library(data.table)
setDTthreads(1)

# record keeping ----------------------------------------------------------

score_model = function(method = "unnamed method", graph, true_graph, individual, 
                       simulation, covariate, p, vb_weight = 1) {
  if (is.null(p)) {
    p = nrow(graph) - 1
  }
  
  est_graph = 1 * (graph > 0.5)
  
  data.table(
    method = method,
    sensitivity = sum(est_graph & true_graph) / sum(true_graph),
    specificity = sum(!est_graph & !true_graph) / sum(!true_graph),
    individual = individual, 
    simulation = simulation,
    covariate = covariate,
    p = p,
    w = w
  )
}

score_model = function(method_name = "unnamed method", )

# Parallel control --------------------------------------------------------
library(doRNG)
registerDoFuture()
plan(multisession, workers = parallel::detectCores() / 2)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

progressr::handlers("progress")

# Discrete covariate, independent ------------------------------------------------------
set.seed(1)
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

# covariate matrix
Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = 1, byrow = FALSE)
tau = .1  # weighting bandwidth

# true graph
true_graph = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    true_graph[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(true_graph) = 0
true_individual_graphs = replicate(n, true_graph, simplify = FALSE)


# hyperparameter grid
p_incl_grid = seq(.05, .8, .05)
v_noise = 1
v_slab = 1

n_sim = 50

cvx_wgts = c(.25, .5, .75)

set.seed(1)
progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  foreach(sim_idx = 1:n_sim) %dorng% {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate the data for this iteration
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
      
      # fit the 2 x p distinct regression models
      varbvs_vopt_result = 
        epwpl::wpl_varbvs(data_mat, Z, tau, 
                          v_noise_grid = v_noise, 
                          v_slab_grid = v_slab,
                          p_incl_grid = p_incl_grid,
                          opt = TRUE)
      
      ep_vopt_result = 
        epwpl::wpl_ep(data_mat, Z, tau,
                      v_noise_grid = v_noise,
                      v_slab_grid = v_slab,
                      p_incl_grid = p_incl_grid,
                      damping = .9, k = .99,
                      opt = TRUE)
      
      epwpl::score_graphs(ep_vopt_result, list(true_graph, true_graph))
      epwpl::score_graphs(varbvs_vopt_result, list(true_graph, true_graph))
      
      metrics = rbind(
        score_model("varbvs_fix", mean_symmetrize(vsvb_result$graphs[[1]]), true_graph,
                    1, sim_idx, -.1, p, 1),
        score_model("varbvs_fix", mean_symmetrize(vsvb_result$graphs[[n]]), true_graph,
                    2, sim_idx, .1, p, 1),
        
        score_model("EP_fix", mean_symmetrize(ep_fix_result$graphs[[1]]), true_graph, 
                    1, sim_idx, -.1, p, 0),
        score_model("EP_fix", mean_symmetrize(ep_fix_result$graphs[[n]]), true_graph, 
                    2, sim_idx, .1, p, 0),
        
        score_model("EP_opt", mean_symmetrize(ep_opt_result$graphs[[1]]), true_graph, 
                    1, sim_idx, -.1, p, 0),
        score_model("EP_opt", mean_symmetrize(ep_opt_result$graphs[[n]]), true_graph, 
                    2, sim_idx, .1, p, 0)
        
        score_model("combo_25", 
                    mean_symmetrize(.25 * vsvb_result$graphs[[1]] + .75 * ep_fix_result$graphs[[1]]))
      )
      
      combo_metrics = rbindlist(lapply(cvx_wgts, function(w) {
        rbind(
          score_model("combo", 
                      mean_symmetrize(w * vsvb_result$graphs[[1]] + 
                                        (1-w) * ep_fix_result$graphs[[1]]),
                      true_graph,
                      1, sim_idx, -.1, p, w),
          score_model("combo",
                      mean_symmetrize(w * vsvb_result$graphs[[n]] + 
                                        (1-w) * ep_fix_result$graphs[[2]]),
                      true_graph,
                      1, sim_idx, .1, p, w)
        )
      }))
      
      # incl_prob = list(
      #   vsvb = vsvb_result,
      #   ep_fix = ep_fix_result,
      #   ep_opt = ep_opt_result
      # )
      # filename = paste0("data/discrete_independent_graph/",
      #                   Sys.Date(),
      #                   "_sim_idx=", sim_idx,
      #                   "_graphs.RDS")
      # saveRDS(incl_prob, file = filename)
      
      filename = paste0("data/discrete_independent/",
                        Sys.Date(), 
                        "_sim_idx" = sim_idx,
                        "_covariate_independent.RDS")
      saveRDS(rbind(metrics, combo_metrics), file = filename)
  }
})
