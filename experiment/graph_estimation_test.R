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

# hyperparameter prior grid
v_noise = 1
v_slab = 3
n_pip = 20

v_noise_grid = rep(v_noise, n_pip)
v_slab_grid = rep(v_slab, n_pip)


# simulate the data for this iteration
X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
data_mat = rbind(X1, X2)

wpl_ep(data_mat, weight_mat)