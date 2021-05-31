# library(future.apply)
library(data.table)
source("experiment/simulation.R")
setDTthreads(1)

# Parallel control --------------------------------------------------------
# library(future.apply)
# library(doFuture)
library(doRng)
registerDoFuture()
plan(multisession, workers = 8)

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
beta = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(beta) = 0

# initial hyperparameter values
sigma0 = 1
p0 = .2
v_slab = 3

n_sim = 50

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  sim_accuracy = rbindlist(
    future_lapply(1:n_sim, function(sim_idx) {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate the data for this iteration
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
      
      # compute weights
      tau = 1  # bandwidth
      D = matrix(1, n, n)
      for(i in 1:n){
        for(j in 1:n){
          D[i, j] = dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
        }
      }
      for(i in 1:n){
        D[, i] = n * (D[, i] / sum(D[, i])) #Scaling the weights so that they add up to n
      }
      
      # fit the n x p regression models
      graphs = wpl_regression(data_mat, D, sigma0, p0, v_slab, n_threads = 1,
                              blas_threads = 1)
      
      # compute accuracy metrics
      metrics = rbindlist(
        foreach(individual = 1:length(graphs)) %do% {
          graph = graphs[[individual]]
          
          # symmetrize estimated graph
          for(i in 1:(p+1)) {
            for(j in i:(p+1)) {
              graph[i, j] = mean(c(graph[i, j], graph[j, i]))
              graph[j, i] = graph[i, j]
            }
          }
          
          est_graph = 1 * (graph > 0.5)
          
          data.table(
            sensitivity = sum(est_graph & beta) / sum(beta),
            specificity = sum(!est_graph & !beta) / sum(!beta),
            individual = individual, 
            simulation = sim_idx
          )
        }
      )
      
      return(metrics)
    }, future.seed = TRUE)
  )
  
})

filename = paste0("data/", Sys.Date(), "_covariate_independent.RDS")
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

D = matrix(1, n, n)  # weights all 1 in no covariate case

# true graph
beta = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(beta) = 0

# initial hyperparameter values
sigma0 = 1
p0 = .2
v_slab = 3

n_sim = 50

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)

  sim_accuracy = rbindlist(
    future_lapply(1:n_sim, function(sim_idx) {
      
      prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
      
      # simulate data
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      data_mat = rbind(X1, X2)
      
      # fit the n x p regression models
      graphs = wpl_regression(data_mat, D, sigma0, p0, v_slab, n_threads = 1,
                              blas_threads = 1)
      
      # compute accuracy metrics
      metrics = rbindlist(
        foreach(individual = 1:length(graphs)) %do% {
          
          graph = graphs[[individual]]
          
          # symmetrize estimated graph
          for(i in 1:(p+1)) {
            for(j in i:(p+1)) {
              graph[i, j] = mean(c(graph[i, j], graph[j, i]))
              graph[j, i] = graph[i, j]
            }
          }
          
          est_graph = 1 * (graph > 0.5)
          
          data.table(
            sensitivity = sum(est_graph & beta) / sum(beta),
            specificity = sum(!est_graph & !beta) / sum(!beta),
            individual = individual, 
            simulation = sim_idx
          )
        }
      )
      
      return(metrics)
    }, future.seed = TRUE)
  )
})

filename = paste0("data/", Sys.Date(), "_no_covariate.RDS")
saveRDS(sim_accuracy, file = filename)


# Discrete covariate, dependent -------------------------------------------
set.seed(1)
n = 100
n_sim = 50

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)
  
  for (p in c(10, 30, 50)) {
    
    Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
    Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5
    
    Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
    Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2
    
    # covariate matrix
    Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)
    
    # true graphs
    beta_neg = matrix(0, p+1, p+1)
    for(i in 1:(p+1)){
      for(j in 1:(p+1)){
        beta_neg[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
      }}
    diag(beta_neg) = 0
    
    beta_pos = matrix(0, nrow = p+1, ncol = p+1)
    for(i in 1:(p+1)){
      for(j in 1:(p+1)){
        beta_pos[i,j] = (Lam2[i] != 0 & Lam2[j] != 0)
      }}
    diag(beta_pos) = 0
    
    # initial hyperparameter values
    sigma0 = 1
    p0 = .2
    v_slab = 3 
    tau = 1  # bandwidth
    
    n_sim = 50
    
    sim_accuracy = rbindlist(
      future_lapply(1:n_sim, function(sim_idx) {
        
        prog(sprintf("Dimension %g, Simulation %g, %s", p, sim_idx, Sys.time()))
        
        X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
        X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
        data_mat = rbind(X1, X2)
        
        # compute weights
        D = matrix(1, n, n)
        for(i in 1:n){
          for(j in 1:n){
            D[i, j] = dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
          }
        }
        for(i in 1:n){
          D[, i] = n * (D[, i] / sum(D[, i])) # Scaling the weights so that they add up to n
        }
        
        # fit the n x p regression models
        graphs = wpl_regression(data_mat, D, sigma0, p0, v_slab, n_threads = 1,
                                blas_threads = 1, woodbury = FALSE)
        
        # compute accuracy metrics
        metrics = rbindlist(
          lapply(graphs, function(graph) {
            
            # symmetrize estimated graph
            for(i in 1:(p+1)) {
              for(j in i:(p+1)) {
                graph[i, j] = max(graph[i, j], graph[j, i])
                graph[j, i] = graph[i, j]
              }
            }
            
            est_graph = 1 * (graph > 0.5)
            
            if (i <= (n/2)) 
              beta = beta_neg
            else
              beta = beta_pos
            
            data.table(
              sensitivity = sum(est_graph & beta) / sum(beta),
              specificity = sum(!est_graph & !beta) / sum(!beta),
              individual = i, 
              simulation = sim_idx,
              p = p
            )
          })
        )
        
        return(metrics)
      }, future.seed = TRUE)
    )
    
    filename = paste0("data/", Sys.Date(), "_p=", p, "_dependent_covariate.RDS")
    saveRDS(sim_accuracy, file = filename)
  }
})

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
      
      graphs = wpl_regression(X, D, sigma0, p0, v_slab, n_threads = 1,
                              blas_threads = 1, woodbury = FALSE)
      
      incl_prob = vector("list", length(graphs))
      
      for (i in 1:length(graphs)) {
        
        graph = graphs[[i]]
        
        for(k in 1:(p+1)) {
          for(j in k:(p+1)) {
            graph[k, j] = mean(c(graph[k, j], graph[j, k]))
            graph[j, k] = graph[k, j]
          }
        }
        
        incl_prob[[i]] = data.table(
          incl_prob_12 = graph[1, 2],
          incl_prob_13 = graph[1, 3],
          individual = i,
          simulation = sim_idx
        )
      }
      
      return(rbindlist(incl_prob))
    }
})

filename = paste0("data/", Sys.Date(), "_continuous_covariate.RDS")
saveRDS(incl_prob_dt, file = filename)

