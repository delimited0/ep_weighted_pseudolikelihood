library(data.table)
source("experiment/simulation.R")
setDTthreads(1)

# record keeping ----------------------------------------------------------

score_model = function(method = "unnamed method", graph, true_graph, individual, 
                       simulation, covariate, p, w = 1) {
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
Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)

# true graph
true_graph = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    true_graph[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(true_graph) = 0

# compute weights
tau = .1  # bandwidth
weight_mat = epwpl::weight_matrix(n, Z, tau)

# only two covariate levels --> only two weightings
weight_mat_fit = weight_mat[c(1, n), ]  

# hyperparameter grid
p_incl_grid = seq(.05, .25, .05)
n_grid = length(p_incl_grid)
# v_noise_grid = rep(1, n_grid)
# v_slab_grid = rep(3, n_grid)

v_noise_grid = c(0.01, 0.05, 0.1, 0.5, 1)
v_slab_grid = c(0.01, 0.05, 0.1, 0.5, 1)

# averaging over cartesian product of hyper param settings
theta = expand.grid(list(
  p_incl_grid = p_incl_grid, 
  v_noise = v_noise_grid,
  v_slab = v_slab_grid
))

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
      vb_result = wpl_varbvs(data_mat, weight_mat_fit, 
                             v_noise_grid = theta$v_noise_grid, 
                             v_slab_grid = theta$v_slab_grid,
                             p_incl_grid = theta$p_incl_grid,
                             opt = FALSE)
      
      ep_result = wpl_ep_gss(data_mat, weight_mat_fit,
                             v_noise_grid,
                             v_slab_grid,
                             p_incl_grid, 
                             damping = 1, k = .99,
                             opt = FALSE,
                             opt_upper = c(100, Inf))
      
      # vb_result = wpl_vb_regression(data_mat, weight_mat_fit, 
      #                               ep_result$sigma_noise, p0, ep_result$v_slab)
      
      metrics = rbind(
        score_model("VB", mean_symmetrize(vsvb_result$graphs[[1]]), true_graph,
                    1, sim_idx, -.1, p, 1),
        score_model("VB", mean_symmetrize(vsvb_result$graphs[[n]]), true_graph,
                    2, sim_idx, .1, p, 1),
        
        score_model("EP_fix", mean_symmetrize(ep_fix_result$graphs[[1]]), true_graph, 
                    1, sim_idx, -.1, p, 0),
        score_model("EP_fix", mean_symmetrize(ep_fix_result$graphs[[2]]), true_graph, 
                    2, sim_idx, .1, p, 0),
        
        score_model("EP_opt", mean_symmetrize(ep_opt_result$graphs[[1]]), true_graph, 
                    1, sim_idx, -.1, p, 0),
        score_model("EP_opt", mean_symmetrize(ep_opt_result$graphs[[2]]), true_graph, 
                    2, sim_idx, .1, p, 0)
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
