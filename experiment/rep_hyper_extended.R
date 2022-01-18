# Command line arguments --------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
output_path = args[1]
p_start = as.numeric(args[2])
p_stop = as.numeric(args[3])
p_by = as.numeric(args[4])
workers = as.numeric(args[5])

print(paste0("Writing results to ", output_path))
print(paste0("Averaging over pip values ", paste0(seq(p_start, p_stop, p_by), collapse=",")))

# Parallel control --------------------------------------------------------
library(doRNG)
library(doFuture)
registerDoFuture()
plan(multisession, workers = workers)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

library(data.table)
setDTthreads(1)

options(progressr.enable = TRUE)
progressr::handlers("progress")

# Discrete covariate, independent -----------------------------------------
n = 100
p = 10
n_nonzero = 3
tau = .1
v_noise = 1
v_slab = 3
n_sim = 50

ep_damping = .95

graph_dims = c(10, 30, 50, 100)
p_incl_grids = lapply(graph_dims, function(d) {
  true_p = n_nonzero / d
  seq(true_p * 1/3, true_p * 5/3, length.out = 11)
})
covariate_dims = c(1, 5)
covariates = list(c(-.1, .1), c(-.2, 0, .2), c(-.3, -.1, .1, .3))

p_incl_grid_map = data.table(graph_dim = graph_dims, p_incl_grid = p_incl_grids)
setkey(p_incl_grid_map, graph_dim)

settings = CJ(
  graph_dim = graph_dims,
  covariate_dim = covariate_dims,
  covariate = covariates,
  sim_idx = 1:n_sim,
  sorted = FALSE
)
setkey(settings, graph_dim)
settings = merge(settings, p_incl_grid_map)

set.seed(1)
metrics = foreach(i = 1:nrow(settings), .combine = rbind, .errorhandling = 'remove') %dorng% {
  prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
  
  setting = settings[i, ]
  p = setting$graph_dim
  cov_dim = setting$covariate_dim
  covariate = setting$covariate[[1]]
  p_incl_grid = setting$p_incl_grid[[1]]
  
  lambda = c(rep(15, 4), rep(0, p-3))
  nugget = 10
  Var = solve(lambda %*% t(lambda) + diag(rep(nugget, p+1)))
  
  n_cov = length(covariate)
  covariates_full = matrix(rep(covariate, each = n / n_cov),
                           nrow = n, ncol = cov_dim, byrow = FALSE)
  covariates_levels = matrix(covariate, nrow = n_cov, ncol = cov_dim, byrow = FALSE)
  weight_idx = seq(1, n * (n_cov-1) / n_cov + 1, n / n_cov)
  weight_mat = epwpl::weight_matrix(covariates_full, tau)[weight_idx, ]
  
  # true graph
  true_graph = matrix(0, p+1, p+1)
  for(i in 1:(p+1)){
    for(j in 1:(p+1)){
      true_graph[i,j] = (lambda[i] != 0 & lambda[j] != 0)
    }}
  diag(true_graph) = 0
  true_individual_graphs = replicate(n_cov, true_graph, simplify = FALSE)
  
  # simulate the data for this iteration
  data_mat = MASS::mvrnorm(n, rep(0, p+1), Var)
  
  # fit the 2 x p distinct regression models
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
                  damping = ep_damping, k = .99,
                  opt = TRUE,
                  woodbury = FALSE,
                  lb = 1e-6,
                  opt_method = "Nelder-Mead")
  
  lapply(ep_vopt_result$individuals, function(indiv) {
    colSums(sapply(indiv$fits, function(fit) {
      fit$mliks * fit$w
    }))
  })
  
  lapply(varbvs_vopt_result$individuals, function(indiv) {
    colSums(sapply(indiv$fits, function(fit) {
      fit$logw * fit$w
    }))
  })
  
  vb_weight = 
 
  metrics = rbind(
    epwpl::score_graphs(ep_vopt_result, true_individual_graphs),
    epwpl::score_graphs(varbvs_vopt_result, true_individual_graphs),
    
    epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .25, true_individual_graphs),
    epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .5, true_individual_graphs),
    epwpl::score_combo_graphs(ep_vopt_result, varbvs_vopt_result, .75, true_individual_graphs)
  )
  metrics[, sim_idx := sim_idx] 
}
