n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

# covariate matrix
Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)

# compute weights
tau = 1
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

X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
data_mat = rbind(X1, X2)

# initial hyperparameter values
sigma0 = 1
p0 = .2
v_slab = 3 

# fit the n x p regression models
graphs = wpl_regression(data_mat, D[c(1, 100), ], sigma0, p0, v_slab, n_threads = 2,
                        blas_threads = 1, woodbury = FALSE)

neg_graph = graphs[[1]]

# symmetrize estimated graph
for(i in 1:(p+1)) {
  for(j in i:(p+1)) {
    neg_graph[i, j] = max(neg_graph[i, j], neg_graph[j, i])
    neg_graph[j, i] = neg_graph[i, j]
  }
}

est_neg_graph = 1 * (neg_graph > 0.5)


