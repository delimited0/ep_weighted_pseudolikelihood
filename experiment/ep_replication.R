# library(future.apply)
library(data.table)
source("experiment/simulation.R")
setDTthreads(1)

# Parallel control --------------------------------------------------------
# library(future.apply)
# library(doFuture)
library(doRNG)
registerDoFuture()
plan(multisession, workers = 8)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

progressr::handlers("progress")

# RhpcBLASctl::blas_set_num_threads(1)

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
tau = 1  # bandwidth
weight_mat = weight_matrix(n, Z)

# only two covariate levels --> only two weightings
weight_mat_fit = weight_mat[c(1, n), ]  

# initial hyperparameter values
sigma0 = 1
p0 = .2
v_slab = 3

n_sim = 50

cvx_wgts = c(.25, .5, .75)

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  sim_accuracy = rbindlist(
    foreach(sim_idx = 1:n_sim) %dorng% {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate the data for this iteration
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
  
      # fit the 2 x p distinct regression models
      
      
      # run sutanoy's vsvb first
      vsvb_result =
        wpl_vsvb_regression(data_mat, weight_mat,
                            sigma0, p0, v_slab, tune = TRUE)
      
      # use the the hyperparams selected there to fit ep
      ep_fix_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                                       sqrt(vsvb_result$sigma0sq), 
                                       p0,
                                       vsvb_result$v_slab, opt = FALSE)
      
      # using Nelder Mead hyperparam optimization
      ep_opt_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                                        sigma0, 
                                        p0,
                                        v_slab, 
                                        opt = TRUE)
      
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
  )
})

filename = paste0("data/discrete_independent/", Sys.Date(), "_covariate_independent.RDS")
saveRDS(sim_accuracy, file = filename)


# No covariate model ------------------------------------------------------
set.seed(1)
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

# covariate matrix
Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)

tau = 1  # bandwidth

weight_mat = matrix(1, n, n)  # weights all 1 in no covariate case
weight_mat_fit = matrix(1, 2, n)

# true graph
true_graph = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    true_graph[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(true_graph) = 0

# initial hyperparameter values
sigma0 = 1
p0 = .2
v_slab = 3

n_sim = 50

cvx_wgts = c(.25, .5, .75)

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)

  sim_accuracy = rbindlist(
    foreach(sim_idx = 1:n_sim) %dorng% {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate data
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
      
      # fit the 2 x p distinct regression models
      
      # run sutanoy's vsvb first
      vsvb_result =
        wpl_vsvb_regression(data_mat, weight_mat,
                            sigma0, p0, v_slab, tune = TRUE)
      
      # use the the hyperparams selected there to fit ep
      ep_fix_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                                        sqrt(vsvb_result$sigma0sq), 
                                        p0,
                                        vsvb_result$v_slab, opt = FALSE)
      
      # using Nelder Mead hyperparam optimization
      ep_opt_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                                        sigma0, 
                                        p0,
                                        v_slab, 
                                        opt = TRUE)
      
      # compute accuracy metrics
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
      
      filename = paste0("data/no_covariate/",
                        Sys.Date(), 
                        "_sim_idx" = sim_idx,
                        "no_covariate.RDS")
      saveRDS(rbind(metrics, combo_metrics), file = filename)
      
      # return(metrics)
    }
  )
})

filename = paste0("data/no_covariate/", Sys.Date(), "_no_covariate.RDS")
saveRDS(sim_accuracy, file = filename)


# Discrete covariate, dependent -------------------------------------------
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
  Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)
  
  # compute weights
  tau = 1  # bandwidth
  weight_mat = weight_matrix(n, Z)
  weight_mat_fit = weight_mat[c(1, n), ]
  
  # true graphs
  true_graph_neg = Prec1 != 0
  diag(true_graph_neg) = 0
  true_graph_pos = Prec2 != 0
  diag(true_graph_pos) = 0
  
  # initial hyperparameter values
  sigma0 = 1
  p0 = .2
  v_slab = 3 
  
  n_sim = 50
  
  progressr::with_progress({
    prog = progressr::progressor(along = 1:n_sim)
    
    sim_accuracy = 
      foreach(sim_idx = 1:n_sim) %dorng% {
        
        prog(sprintf("Dimension %g, Simulation %g, %s", p, sim_idx, Sys.time()))
        
        X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
        X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
        data_mat = rbind(X1, X2)
        
        # fit the n x p regression models
        
        # run sutanoy's vsvb first
        vsvb_result =
          wpl_vsvb_regression(data_mat, weight_mat,
                              sigma0, p0, v_slab, tune = TRUE)
        
        # use the the hyperparams selected there to fit ep
        ep_fix_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                                          sqrt(vsvb_result$sigma0sq), 
                                          p0,
                                          vsvb_result$v_slab, opt = FALSE)
        
        # using Nelder Mead hyperparam optimization
        ep_opt_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                                          sigma0, 
                                          p0,
                                          v_slab, 
                                          opt = TRUE)
        
        # compute accuracy metrics
        metrics = rbind(
          score_model("VB", mean_symmetrize(vsvb_result$graphs[[1]]), true_graph_neg,
                      1, sim_idx, -.1, p),
          score_model("VB", mean_symmetrize(vsvb_result$graphs[[n]]), true_graph_pos,
                      2, sim_idx, .1, p),
          
          score_model("EP_fix", mean_symmetrize(ep_fix_result$graphs[[1]]), true_graph_neg, 
                      1, sim_idx, -.1, p),
          score_model("EP_fix", mean_symmetrize(ep_fix_result$graphs[[2]]), true_graph_pos, 
                      2, sim_idx, .1, p),
          
          score_model("EP_opt", mean_symmetrize(ep_opt_result$graphs[[1]]), true_graph_neg, 
                      1, sim_idx, -.1, p),
          score_model("EP_opt", mean_symmetrize(ep_opt_result$graphs[[2]]), true_graph_pos, 
                      2, sim_idx, .1, p)
        )
        
        filename = paste0("data/discrete_dependent/",
                          Sys.Date(),
                          "_sim=", sim_idx,
                          "_p=", p, "_dependent_covariate.RDS")
        saveRDS(metrics, file = filename)
        
        return(metrics)
      }
  })
}

# Continuous covariate ----------------------------------------------------

n = 180
p = 4
MAXITER = 1
STR = 1
in_pr_13 = matrix(0, MAXITER, n)
in_pr_12 = in_pr_13

Var_cont = function(z) {
  
  pr = matrix(0, p+1, p+1)
  diag(pr) = 2
  
  pr[2,3] = STR
  pr[1,2] = STR*((z>-1) && (z< -.33)) + (STR - STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (0)*((z>0.43) && (z<1))
  pr[1,3] = 0*((z>-1) && (z< -.33)) + (STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (STR)*((z>0.43) && (z<1))
  
  pr[2,1] = pr[1,2]
  pr[3,1] = pr[1,3]
  pr[3,2] = pr[2,3]
  
  
  Var = solve(pr)
  return(Var)
}

sensitivity_20 = matrix(0, MAXITER, 1)
specificity_20 = sensitivity_20
sensitivity_90 = sensitivity_20
specificity_90 = sensitivity_20
sensitivity_160 = sensitivity_20
specificity_160 = sensitivity_20

Z = c(seq(-0.99, -0.331, (-.331+.99)/59), 
      seq(-0.229,0.329,(.329+.229)/59),
      seq(0.431,.99,(.99-.431)/59))
Z = matrix(Z, n, 1)
X = matrix(0, n, p+1)

p0 = .2
v_slab = 3
sigma0 = 1
tau = 0.56

n_sim = 50

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  # incl_prob_dt = rbindlist(
    
  incl_prob_dt = 
    foreach(sim_idx = 1:n_sim, .combine = rbind) %dorng% {
      
      prog(sprintf("Simulation %g, %s", p, sim_idx, Sys.time()))
      
      # simulate data
      for(i in 1:n) {
        X[i, ] = MASS::mvrnorm(1, rep(0, p+1), Var_cont(Z[i]))
      }
      data_mat = X
      
      # weight matrix
      
      D = matrix(1, n, n)
      for(i in 1:n) {
        for(j in 1:n) {
          D[j, i] = dnorm(norm(Z[i,]-Z[j,],"2"), 0, tau)
        }
      }
      for(i in 1:n){
        D[, i] = n*(D[,i] / sum(D[,i]))
      }
      weight_mat = D
      
      # run sutanoy's vsvb first
      vsvb_result =
        wpl_vsvb_regression(data_mat, weight_mat,
                            sigma0, p0, v_slab, tune = TRUE)
      
      # use the the hyperparams selected there to fit ep
      ep_fix_result = wpl_ep_regression(data_mat, weight_mat, 
                                        sqrt(vsvb_result$sigma0sq), 
                                        p0,
                                        vsvb_result$v_slab, opt = FALSE)
      
      # using Nelder Mead hyperparam optimization
      ep_opt_result = wpl_ep_regression(data_mat, weight_mat, 
                                        sigma0, 
                                        p0,
                                        v_slab, 
                                        opt = TRUE)
      
      
      vsvb_probs = rbindlist(lapply(1:n, function(indiv) {
        graph = vsvb_result$graphs[[indiv]]
        incl_prob_model("VB", mean_symmetrize(graph), indiv, sim_idx, Z[indiv, ], p) 
      }))
      
      ep_fix_probs = rbindlist(lapply(1:n, function(i) {
        graph = ep_fix_result$graphs[[i]]
        incl_prob_model("EP_fix", mean_symmetrize(graph), i, sim_idx, Z[i, ], p) 
      }))
      
      ep_opt_probs = rbindlist(lapply(1:n, function(i) {
        graph = ep_opt_result$graphs[[i]]
        incl_prob_model("EP_opt", mean_symmetrize(graph), i, sim_idx, Z[i, ], p) 
      }))

      probs = rbind(vsvb_probs, ep_fix_probs, ep_opt_probs)
      
      filename = paste0("data/continuous/",
                        Sys.Date(),
                        "_sim=", sim_idx,
                        "_p=", p, "_dependent_covariate.RDS")
      saveRDS(probs, file = filename)
    }
})

filename = paste0("data/", Sys.Date(), "_continuous_covariate.RDS")
saveRDS(incl_prob_dt, file = filename)

