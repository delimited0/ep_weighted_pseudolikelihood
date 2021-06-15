library(data.table)
source("experiment/simulation.R")
setDTthreads(1)

library(doRNG)
registerDoFuture()
plan(multisession, workers = 8)

progressr::handlers("progress")

library(ggplot2)

# More distinct covariates ------------------------------------------------

set.seed(1)
n = 100
n_sim = 50

pos_covariate = 1
neg_covariate = -1


progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  for (p in c(10, 30, 50)) {
  # for (p in c(50)) {
    
    Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
    Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5
    
    Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
    Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2
    
    # covariate matrix
    Z = matrix(neg_covariate*(1:n <= n/2) +
               pos_covariate*(1:n > n/2), 
               nrow = n, ncol = p, byrow = FALSE)
    
    # compute weights
    tau = 1  # bandwidth
    weight_mat = weight_matrix(n, Z)
    weight_mat_fit = weight_mat[c(1, n), ]
    
    # true graphs
    true_graph_neg = matrix(0, p+1, p+1)
    for(i in 1:(p+1)){
      for(j in 1:(p+1)){
        true_graph_neg[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
      }}
    diag(true_graph_neg) = 0
    
    true_graph_pos = matrix(0, nrow = p+1, ncol = p+1)
    for(i in 1:(p+1)){
      for(j in 1:(p+1)){
        true_graph_pos[i,j] = (Lam2[i] != 0 & Lam2[j] != 0)
      }}
    diag(true_graph_pos) = 0
    
    # initial hyperparameter values
    sigma0 = 1
    p0 = .2
    v_slab = 3 
    
    n_sim = 50
    
    all_metrics = 
      # lapply(1:n_sim, function(sim_idx) {
      # foreach(sim_idx = 1:n_sim) %do% {
      foreach(sim_idx = 1:n_sim) %dorng% {
        
        prog(sprintf("Dimension %g, Simulation %g, %s", p, sim_idx, Sys.time()))
        
        X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
        X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
        data_mat = rbind(X1, X2)
        
        # fit the n x p regression models
        graphs = wpl_regression(data_mat, weight_mat_fit, sigma0, p0, v_slab, n_threads = 1,
                                blas_threads = 1, woodbury = FALSE)
        
        # compute accuracy metrics
        
        metrics = rbind(
          score_model(mean_symmetrize(graphs[[1]]), true_graph_neg, 1, sim_idx, neg_covariate, p),
          score_model(mean_symmetrize(graphs[[2]]), true_graph_pos, 2, sim_idx, pos_covariate, p)
        )
        
        filename = paste0("data/discrete_dependent_checking/", 
                          Sys.Date(), "_covrange=", pos_covariate, 
                          "_sim=", sim_idx,
                          "_p=", p, "_dependent_covariate.RDS")
        saveRDS(metrics, file = filename)
        
        return(metrics)
      # })
      }
  }
})

dat_names = dir("data/discrete_dependent_checking/wide_range/", full.names = TRUE)

wide_sim = rbindlist(lapply(dat_names, function(filename) readRDS(filename)))
wide_sim_long = melt(wide_sim, id.vars = c("individual", "simulation", "covariate", "p"))

plt = ggplot(wide_sim_long, aes(x = value, fill = as.factor(covariate))) +
  geom_boxplot() +
  facet_grid(rows = vars(variable), cols = vars(p)) +
  coord_flip() + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "Covariate")

# 3 covariate levels ------------------------------------------------------

set.seed(1)
n = 102
n_sim = 50

covariates = c(-.1, 0, .1)

for (p in c(10, 30, 50)) {
  
  Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
  Lam2 = c(rep(0, 4), 3, 3, 3, 3, rep(0, p - 7)) * 5
  Lam3 = c(rep(0, p-3), 3, 3, 3, 3) * 5
  
  Prec1 = Lam1 %*% t(Lam1) + diag(rep(10, p+1))
  Prec2 = Lam2 %*% t(Lam2) + diag(rep(10, p+1))
  Prec3 = Lam3 %*% t(Lam3) + diag(rep(10, p+1))
  
  Var1 = solve(Prec1) 
  Var2 = solve(Prec2) 
  Var3 = solve(Prec3)
  
  # covariate matrix
  Z = matrix(
    covariates[1] * (1:n <= n/3) + covariates[3] * (1:n > (2*n)/3), 
    nrow = n, ncol = p, byrow = FALSE)
  
  # compute weights
  tau = 1  # bandwidth
  weight_mat = weight_matrix(n, Z)
  weight_mat_fit = weight_mat[c(1, (n/2), n), ]
  
  # true graphs
  tg1 = Prec1 != 0
  diag(tg1) = 0
  tg2 = Prec2 != 0
  diag(tg2) = 0
  tg3 = Prec3 != 0
  diag(tg3) = 0
  
  # initial hyperparameter values
  sigma0 = 1
  p0 = .2
  v_slab = 3 
  
  n_sim = 50
  
  progressr::with_progress({
    prog = progressr::progressor(along = 1:n_sim)
    
    all_metrics = 
      # lapply(1:n_sim, function(sim_idx) {
      # foreach(sim_idx = 1:n_sim) %do% {
      foreach(sim_idx = 1:n_sim) %dorng% {
        
        prog(sprintf("Dimension %g, Simulation %g, %s", p, sim_idx, Sys.time()))
        
        X1 = MASS::mvrnorm(n/3, rep(0, p+1), Var1)
        X2 = MASS::mvrnorm(n/3, rep(0, p+1), Var2)
        X3 = MASS::mvrnorm(n/3, rep(0, p+1), Var3)
        data_mat = rbind(X1, X2, X3)
        
        # fit the n x p regression models
        graphs = wpl_regression(data_mat, weight_mat_fit, sigma0, p0, v_slab, n_threads = 1,
                                blas_threads = 1, woodbury = FALSE)
        
        # compute accuracy metrics
        
        metrics = rbind(
          score_model(mean_symmetrize(graphs[[1]]), tg1, 1, sim_idx, covariates[1], p),
          score_model(mean_symmetrize(graphs[[2]]), tg2, 2, sim_idx, covariates[2], p),
          score_model(mean_symmetrize(graphs[[3]]), tg3, 3, sim_idx, covariates[3], p)
        )
        
        filename = paste0("data/discrete_dependent_checking/three_levels/", 
                          Sys.Date(), 
                          "_sim=", sim_idx,
                          "_p=", p, "_dependent_covariate.RDS")
        saveRDS(metrics, file = filename)
        
        return(metrics)
        # })
      }
  })
}

dat_names = dir("data/discrete_dependent_checking/three_levels/", full.names = TRUE)

three_level_sim = rbindlist(lapply(dat_names, function(filename) readRDS(filename)))
three_level_sim_long = melt(three_level_sim, id.vars = c("individual", "simulation", "covariate", "p"))

plt = ggplot(three_level_sim_long, aes(x = value, fill = as.factor(covariate))) +
  geom_boxplot() +
  facet_grid(rows = vars(variable), cols = vars(p)) +
  coord_flip() + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "Covariate")


# 4 covariate levels ------------------------------------------------------

set.seed(1)
n = 100
n_sim = 50

covariates = c(-.5, -.25, .25, .5)


for (p in c(16, 30, 50)) {
  
  Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
  Lam2 = c(rep(0, 4), 3, 3, 3, 3, rep(0, p - 7)) * 5
  Lam3 = c(rep(0, 8), 3, 3, 3, 3, rep(0, p - 11)) * 5
  Lam4 = c(rep(0, p-3), 3, 3, 3, 3) * 5
  
  Prec1 = Lam1 %*% t(Lam1) + diag(rep(10, p+1))
  Prec2 = Lam2 %*% t(Lam2) + diag(rep(10, p+1))
  Prec3 = Lam3 %*% t(Lam3) + diag(rep(10, p+1))
  Prec4 = Lam4 %*% t(Lam4) + diag(rep(10, p+1))
  
  Var1 = solve(Prec1) 
  Var2 = solve(Prec2)
  Var3 = solve(Prec3)
  Var4 = solve(Prec4)
  
  # covariate matrix
  Z = matrix(
    covariates[1] * (1:n <= n/4) + 
      covariates[2] * (1:n > n/4 & 1:n <= n/2) +
      covariates[3] * (1:n > (n/2) & 1:n <= (3*n/4)) +
      covariates[4] * (1:n > (3*n/4)),                         
    nrow = n, ncol = p, byrow = FALSE)
  
  # compute weights
  tau = 1  # bandwidth
  weight_mat = weight_matrix(n, Z)
  weight_mat_fit = weight_mat[c(1, (n/2 + 1), (3*n/4 + 1), n), ]
  
  # true graphs
  tg1 = Prec1 != 0
  diag(tg1) = 0
  tg2 = Prec2 != 0
  diag(tg2) = 0
  tg3 = Prec3 != 0
  diag(tg3) = 0
  tg4 = Prec4 != 0
  diag(tg4) = 0
  
  # initial hyperparameter values
  sigma0 = 1
  p0 = .2
  v_slab = 3 
  
  n_sim = 50
  
  progressr::with_progress({
    prog = progressr::progressor(along = 1:n_sim)
    
    all_metrics = 
      # lapply(1:n_sim, function(sim_idx) {
      # foreach(sim_idx = 1:n_sim) %do% {
      foreach(sim_idx = 1:n_sim) %dorng% {
        
        prog(sprintf("Dimension %g, Simulation %g, %s", p, sim_idx, Sys.time()))
        
        X1 = MASS::mvrnorm(n/4, rep(0, p+1), Var1)
        X2 = MASS::mvrnorm(n/4, rep(0, p+1), Var2)
        X3 = MASS::mvrnorm(n/4, rep(0, p+1), Var3)
        X4 = MASS::mvrnorm(n/4, rep(0, p+1), Var4)
        data_mat = rbind(X1, X2, X3, X4)
        
        # fit the n x p regression models
        graphs = wpl_regression(data_mat, weight_mat_fit, sigma0, p0, v_slab, n_threads = 1,
                                blas_threads = 1, woodbury = FALSE)
        
        # compute accuracy metrics
        metrics = rbind(
          score_model(mean_symmetrize(graphs[[1]]), tg1, 1, sim_idx, covariates[1], p),
          score_model(mean_symmetrize(graphs[[2]]), tg2, 2, sim_idx, covariates[2], p),
          score_model(mean_symmetrize(graphs[[3]]), tg3, 3, sim_idx, covariates[3], p),
          score_model(mean_symmetrize(graphs[[4]]), tg4, 4, sim_idx, covariates[4], p)
        )
        
        filename = paste0("data/discrete_dependent_checking/four_levels/", 
                          Sys.Date(), 
                          "_sim=", sim_idx,
                          "_p=", p, "_dependent_covariate.RDS")
        saveRDS(metrics, file = filename)
        
        return(metrics)
        # })
      }
  })
}

dat_names = dir("data/discrete_dependent_checking/four_levels/", full.names = TRUE)

four_level_sim = rbindlist(lapply(dat_names, function(filename) readRDS(filename)))
four_level_sim_long = melt(four_level_sim, id.vars = c("individual", "simulation", "covariate", "p"))

plt = ggplot(four_level_sim_long, aes(x = value, fill = as.factor(covariate))) +
  geom_boxplot() +
  facet_grid(rows = vars(variable), cols = vars(p)) +
  coord_flip() + 
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(fill = "Covariate")
