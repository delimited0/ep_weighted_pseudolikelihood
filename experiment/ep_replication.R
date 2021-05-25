library(future.apply)
library(data.table)

# Parallel control --------------------------------------------------------

# plan(multisession, workers = 4)

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

sim_accuracy = rbindlist(lapply(1:n_sim, function(sim_idx) {
  
  # fit the n x p regression models
  regressions = lapply(1:(p+1), function(resp_index) {
    
    print(paste0(" --- Covariate ", resp_index, " --- "))
    
    y = data_mat[, resp_index]  
    X = data_mat[, -resp_index]
    
    # select an individual
    incl_prob = sapply(1:n, function(i) {
      
      print(paste0("Individual ", i))
      
      sq_weights = sqrt(D[i, ])
      y_weighted = y * sq_weights
      X_weighted = X * sq_weights
      
      # fit model
      fit = epwpl::ep_wlr(X_weighted, y_weighted, sigma0, p0, v_slab)
      t(plogis(fit$p))
    })
    
    return(incl_prob)
  })
  
  # compute accuracy metrics
  accuracy = rbindlist(
    lapply(1:n, function(i) {
      
      graph_prob = matrix(0, p+1, p+1)
      
      for(j in 1:(p+1)) {
        graph_prob[j, -j] = regressions[[j]][, i] #Individual specific inclusion probability matrix
      }
      
      for(i in 1:(p+1)) {
        for(j in i:(p+1)) {
          graph_prob[i, j] = mean(c(graph_prob[i, j], graph_prob[j, i]))
          graph_prob[j, i] = graph_prob[i, j]
        }
      }
      
      est_graph = 1 * (graph_prob> 0.5)
      
      data.table(
        sensitivity = sum(est_graph & beta) / sum(beta),
        specificity = sum(!est_graph & !beta) / sum(!beta),
        individual = i, 
        simulation = sim_idx
      )
    })
  )
  
  return(accuracy)
}))

filename = paste0("../data/", Sys.Date(), "_covariate_independent.RDS")
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

sim_accuracy = rbindlist(lapply(1:n_sim, function(sim_idx) {
  
  # fit the n x p regression models
  regressions = lapply(1:(p+1), function(resp_index) {
    
    print(paste0(" --- Covariate ", resp_index, " --- "))
    
    y = data_mat[, resp_index]  
    X = data_mat[, -resp_index]
    
    # select an individual
    incl_prob = sapply(1:n, function(i) {
      
      print(paste0("Individual ", i))
      
      sq_weights = sqrt(D[i, ])
      y_weighted = y * sq_weights
      X_weighted = X * sq_weights
      
      # fit model
      fit = epwpl::ep_wlr(X_weighted, y_weighted, sigma0, p0, v_slab)
      t(plogis(fit$p))
    })
    
    return(incl_prob)
  })
  
  # compute accuracy metrics
  accuracy = rbindlist(
    lapply(1:n, function(i) {
      
      graph_prob = matrix(0, p+1, p+1)
      
      for(j in 1:(p+1)) {
        graph_prob[j, -j] = regressions[[j]][, i] #Individual specific inclusion probability matrix
      }
      
      for(i in 1:(p+1)) {
        for(j in i:(p+1)) {
          graph_prob[i, j] = mean(c(graph_prob[i, j], graph_prob[j, i]))
          graph_prob[j, i] = graph_prob[i, j]
        }
      }
      
      est_graph = 1 * (graph_prob> 0.5)
      
      data.table(
        sensitivity = sum(est_graph & beta) / sum(beta),
        specificity = sum(!est_graph & !beta) / sum(!beta),
        individual = i, 
        simulation = sim_idx
      )
    })
  )
  
  return(accuracy)
}))

filename = paste0("../data/", Sys.Date(), "_no_covariate.RDS")
saveRDS(sim_accuracy, file = filename)


# Discrete covariate, dependent -------------------------------------------
set.seed(1)
n = 100
p = 10

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
  D[, i] = n * (D[, i] / sum(D[, i])) #Scaling the weights so that they add up to n
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

sim_accuracy = rbindlist(lapply(1:n_sim, function(sim_idx) {
  
  # fit the n x p regression models
  regressions = lapply(1:(p+1), function(resp_index) {
    
    print(paste0(" --- Covariate ", resp_index, " --- "))
    
    y = data_mat[, resp_index]  
    X = data_mat[, -resp_index]
    
    # select an individual
    incl_prob = sapply(1:n, function(i) {
      
      print(paste0("Individual ", i))
      
      sq_weights = sqrt(D[i, ])
      y_weighted = y * sq_weights
      X_weighted = X * sq_weights
      
      # fit model
      fit = epwpl::ep_wlr(X_weighted, y_weighted, sigma0, p0, v_slab)
      t(plogis(fit$p))
    })
    
    return(incl_prob)
  })
  
  # compute accuracy metrics
  accuracy = rbindlist(
    lapply(1:n, function(i) {
      
      graph_prob = matrix(0, p+1, p+1)
      
      for(j in 1:(p+1)) {
        graph_prob[j, -j] = regressions[[j]][, i] #Individual specific inclusion probability matrix
      }
      
      for(i in 1:(p+1)) {
        for(j in i:(p+1)) {
          graph_prob[i, j] = mean(c(graph_prob[i, j], graph_prob[j, i]))
          graph_prob[j, i] = graph_prob[i, j]
        }
      }
      
      est_graph = 1 * (graph_prob> 0.5)
      
      data.table(
        sensitivity = sum(est_graph & beta) / sum(beta),
        specificity = sum(!est_graph & !beta) / sum(!beta),
        individual = i, 
        simulation = sim_idx
      )
    })
  )
  
  return(accuracy)
}))

filename = paste0("../data/", Sys.Date(), "_dependent_covariate.RDS")
saveRDS(sim_accuracy, file = filename)



