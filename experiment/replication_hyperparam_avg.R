library(data.table)
setDTthreads(1)

# Parallel control --------------------------------------------------------
library(doRNG)
library(doFuture)
registerDoFuture()
plan(multisession, workers = 8)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

progressr::handlers("progress")


# Optimize variance -------------------------------------------------------


## Discrete covariate, independent ------------------------------------------------------
set.seed(1)
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

# covariate matrix
Z_all = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = 1, byrow = FALSE)
Z_2 = matrix(c(-.1, .1), nrow = 2, ncol = 1, byrow = FALSE)
tau = .1  # weighting bandwidth
weight_mat = epwpl::weight_matrix(Z_all, tau)[c(1, n), ]

# true graph
true_graph = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    true_graph[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(true_graph) = 0
true_individual_graphs = replicate(2, true_graph, simplify = FALSE)

# hyperparameter grid
p_incl_grid = seq(.05, .8, .05)
v_noise = 1
v_slab = 3

n_sim = 50

set.seed(1)
progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  metrics = foreach(sim_idx = 1:n_sim, .combine = rbind) %dorng% {
  # for (sim_idx in 1:n_sim) {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate the data for this iteration
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
      
      # fit the 2 x p distinct regression models
      varbvs_vopt_result = 
        epwpl::wpl_varbvs(data_mat, 
                          covariates = Z_2, tau = tau, weight_mat = weight_mat, 
                          v_noise_grid = v_noise, 
                          v_slab_grid = v_slab,
                          p_incl_grid = p_incl_grid,
                          opt = TRUE)
      
      ep_vopt_result = 
        epwpl::wpl_ep(data_mat, 
                      covariates = Z_2, tau = tau, weight_mat = weight_mat,
                      v_noise_grid = v_noise,
                      v_slab_grid = v_slab,
                      p_incl_grid = p_incl_grid,
                      damping = .9, k = .99,
                      opt = TRUE, 
                      opt_method = "Nelder-Mead")
      
      metrics = rbind(
        epwpl::score_graphs(ep_vopt_result, true_individual_graphs),
        epwpl::score_graphs(varbvs_vopt_result, true_individual_graphs),
        epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .25, true_individual_graphs),
        epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .5, true_individual_graphs),
        epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .75, true_individual_graphs)
      )
      metrics[, sim_idx := sim_idx]
      
      metrics
  }
})

filename = paste0("data/hyper_avg/discrete_independent.RDS")
saveRDS(metrics, file = filename)

## No covariate model ------------------------------------------------------
set.seed(1)
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

Z = matrix(NA, 1, 1)
weight_mat = matrix(1, 1, n)  # weights all 1 in no covariate case
tau = .1

# true graph
true_graph = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    true_graph[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(true_graph) = 0
true_individual_graphs = list(true_graph)

# hyperparameter grid
p_incl_grid = seq(.05, .8, .05)
v_noise = 1
v_slab = 3

n_sim = 50

set.seed(1)
progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  metrics = foreach(sim_idx = 1:n_sim, .combine = rbind) %dorng% {
    
    prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
    
    # simulate data
    X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
    X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
    data_mat = rbind(X1, X2)
    
    # fit the 2 x p distinct regression models
    varbvs_vopt_result = 
      epwpl::wpl_varbvs(data_mat, 
                        covariates = Z, tau = tau, weight_mat = weight_mat, 
                        v_noise_grid = v_noise, 
                        v_slab_grid = v_slab,
                        p_incl_grid = p_incl_grid,
                        opt = TRUE)
    
    ep_vopt_result = 
      epwpl::wpl_ep(data_mat, 
                    covariates = Z, tau = tau, weight_mat = weight_mat,
                    v_noise_grid = v_noise,
                    v_slab_grid = v_slab,
                    p_incl_grid = p_incl_grid,
                    damping = .9, k = .99,
                    opt = TRUE, 
                    opt_method = "Nelder-Mead")
    
    metrics = rbind(
      epwpl::score_graphs(ep_vopt_result, true_individual_graphs),
      epwpl::score_graphs(varbvs_vopt_result, true_individual_graphs),
      epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .25, true_individual_graphs),
      epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .5, true_individual_graphs),
      epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .75, true_individual_graphs)
    )
    metrics[, sim_idx := sim_idx]
    
    metrics
  }
})

filename = paste0("data/hyper_avg/no_covariate.RDS")
saveRDS(metrics, file = filename)

## Discrete covariate, dependent ------------------------------------------------------

set.seed(1)
n = 100
n_sim = 50

for (p in c(10, 30, 50)) {
  
  Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
  Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5
  
  Prec1 = Lam1 %*% t(Lam1) + diag(rep(10, p+1))
  Prec2 = Lam2 %*% t(Lam2) + diag(rep(10, p+1))
  
  Var1 = solve(Prec1)
  Var2 = solve(Prec2)
  
  # covariate matrix
  Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = 1, byrow = FALSE)
  Z_all = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)
  Z_2 = matrix(c(-.1, .1), nrow = 2, ncol = 1, byrow = FALSE)
  tau = .1  # weighting bandwidth
  weight_mat = epwpl::weight_matrix(Z_all, tau)[c(1, n), ]
  
  # true graphs
  true_graph_neg = Prec1 != 0
  diag(true_graph_neg) = 0
  true_graph_pos = Prec2 != 0
  diag(true_graph_pos) = 0
  true_individual_graphs = list(true_graph_neg, true_graph_pos)
  
  # hyperparameter grid
  p_incl_grid = seq(.05, .8, .05)
  sigma0 = 1
  v_slab = 3 
  
  set.seed(1)
  progressr::with_progress({
    prog = progressr::progressor(along = 1:n_sim)
    
    metrics = foreach(sim_idx = 1:n_sim, .combine = rbind) %dorng% {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate data
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
      
      # fit the 2 x p distinct regression models
      varbvs_vopt_result = 
        epwpl::wpl_varbvs(data_mat, 
                          covariates = Z_2, tau = tau, weight_mat = weight_mat, 
                          v_noise_grid = v_noise, 
                          v_slab_grid = v_slab,
                          p_incl_grid = p_incl_grid,
                          opt = TRUE)
      
      ep_vopt_result = 
        epwpl::wpl_ep(data_mat, 
                      covariates = Z_2, tau = tau, weight_mat = weight_mat,
                      v_noise_grid = v_noise,
                      v_slab_grid = v_slab,
                      p_incl_grid = p_incl_grid,
                      damping = .9, k = .99,
                      opt = TRUE, 
                      opt_method = "Nelder-Mead")
      
      metrics = rbind(
        epwpl::score_graphs(ep_vopt_result, true_individual_graphs),
        epwpl::score_graphs(varbvs_vopt_result, true_individual_graphs),
        epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .25, true_individual_graphs),
        epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .5, true_individual_graphs),
        epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .75, true_individual_graphs)
      )
      metrics[, sim_idx := sim_idx]
      
      metrics
    }
  })
  
  filename = paste0("data/hyper_avg/discrete_dependent_", "p=", p, ".RDS")
  saveRDS(metrics, file = filename)
}



  

