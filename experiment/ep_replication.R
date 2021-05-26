library(data.table)
setDTthreads(1)

# Parallel control --------------------------------------------------------
# library(future.apply)
library(doFuture)
registerDoFuture()
plan(multisession, workers = 4)
# RhpcBLASctl::blas_set_num_threads(1)

library(progressr)
handlers(global = TRUE)
handlers("progress")

# Discrete covariate, independent ------------------------------------------------------
set.seed(1)
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)

data_mat = rbind(X1, X2)

# covariate matrix
Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)

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
  # D[,i]=1 # When there is no covariate information, set the weights to be 1 throughout.
}

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
      
      # print(paste0(" --- Simulation ", sim_idx, " ---"))
      
      # fit the n x p regression models
      graphs = wpl_regression(data_mat, D, sigma0, p0, v_slab, n_threads = 1,
                              blas_threads = 1)
      
      # compute accuracy metrics
      metrics = rbindlist(
        lapply(graphs, function(graph) {
          
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
            individual = i, 
            simulation = sim_idx
          )
        })
      )
      
      return(metrics)
    })
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

X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)

data_mat = rbind(X1, X2)

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

  sim_accuracy = rbindlist(future_lapply(1:n_sim, function(sim_idx) {
    
    prog(sprintf("Simulation %g, %s", sim_idx, Sys.time()))
    
    # fit the n x p regression models
    graphs = wpl_regression(data_mat, D, sigma0, p0, v_slab, n_threads = 1,
                            blas_threads = 1)
    
    # compute accuracy metrics
    metrics = rbindlist(
      lapply(graphs, function(graph) {
        
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
          individual = i, 
          simulation = sim_idx
        )
      })
    )
    
    return(metrics)
  }))
})

filename = paste0("data/", Sys.Date(), "_no_covariate.RDS")
saveRDS(sim_accuracy, file = filename)


# Discrete covariate, dependent -------------------------------------------
set.seed(1)
n = 100

progressr::with_progress({
  prog = progressr::progressor(along = 1:n_sim)

  disc_cov_accuracy = rbindlist(
    lapply(c(10, 30, 50), function(p) {
      
      Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
      Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5
      
      Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
      Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2
      
      X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
      X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
      
      data_mat = rbind(X1, X2)
      
      # covariate matrix
      Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)
      
      # compute weights
      tau = 1  # bandwidth
      D = matrix(1, n, n)
      for(i in 1:n){
        for(j in 1:n){
          D[i, j] = dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
        }
      }
      for(i in 1:n){
        D[, i] = n * (D[, i] / sum(D[, i])) # Scaling the weights so that they add up to n
      }
      
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
      
      n_sim = 50
      
      sim_accuracy = rbindlist(future_lapply(1:n_sim, function(sim_idx) {
        
        prog(sprintf("Dimension %g, Simulation %g, %s", p, sim_idx, Sys.time()))
        
        # fit the n x p regression models
        graphs = wpl_regression(data_mat, D, sigma0, p0, v_slab, n_threads = 1,
                                blas_threads = 1)
        
        # compute accuracy metrics
        metrics = rbindlist(
          lapply(graphs, function(graph) {
            
            # symmetrize estimated graph
            for(i in 1:(p+1)) {
              for(j in i:(p+1)) {
                graph[i, j] = mean(c(graph[i, j], graph[j, i]))
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
      }))
      
      return(sim_accuracy)
    })
  )
})

filename = paste0("data/", Sys.Date(), "_dependent_covariate.RDS")
saveRDS(sim_accuracy, file = filename)



