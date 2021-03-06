# discrete independent ----

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

# hyperparameter prior grid
v_noise = 1
v_slab = 3
# p_incl_grid = seq(1, p) / (p+1)
# p_incl_grid = seq(.3, .7, .1)
p_incl_grid = seq(.1, .5, .05)
# p_incl_grid = .5
n_pip = length(p_incl_grid)

# v_noise_grid = rep(v_noise, n_pip)
# v_slab_grid = rep(v_slab, n_pip)

# v_noise_grid = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)
# v_slab_grid = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)

v_noise_grid = c(0.01, 0.05, 0.1, 0.5, 1)
v_slab_grid = c(0.01, 0.05, 0.1, 0.5, 1)

# simulate the data for this iteration
X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
data_mat = rbind(X1, X2)


# averaging over cartesian product of hyper param settings
p_incl_grid = seq(.05, .25, .05)
n_pip = length(p_incl_grid)
v_noise_grid = c(0.01, 0.05, 0.1, 0.5, 1)
v_slab_grid = c(0.01, 0.05, 0.1, 0.5, 1)


theta = expand.grid(list(
  p0 = p_incl_grid, 
  v_noise = v_noise_grid,
  v_slab = v_slab_grid
))

tictoc::tic()
result_fix = epwpl::epvbs(data_mat, weight_mat_fit, 
                          theta$v_noise, 
                          theta$v_slab, 
                          theta$p0, opt=FALSE, verbose=FALSE,
                          damping = .9, k = .99,
                        )
tictoc::toc()

tictoc::tic()
result_vb_fix = epwpl::wpl_vb(data_mat, weight_mat_fit,
                              theta$v_noise, 
                              theta$v_slab, 
                              theta$p0, opt=FALSE, verbose=TRUE)
tictoc::tic()

tictoc::tic()
result_varbvs_fix = epwpl::wpl_varbvs(data_mat, weight_mat_fit,
                              theta$v_noise, 
                              theta$v_slab, 
                              theta$p0, opt=FALSE, verbose=TRUE)
tictoc::tic()


result_fix = epwpl::wpl_ep(data_mat, weight_mat_fit, 
                           theta$v_noise, 
                           theta$v_slab, 
                           theta$p0, 
                           damping = .9, k = .99,
                           opt=FALSE, verbose=TRUE)
tictoc::toc()

# optimize over variances, average over only inclusion probabilities
result_opt = epwpl::wpl_ep_gss(data_mat, weight_mat_fit, 
                           v_noise_grid, 
                           v_slab_grid, 
                           p_incl_grid,
                           damping = 1, k = .99,
                           opt=TRUE, verbose=TRUE,
                           opt_upper = c(1000, Inf))

tictoc::tic()
result_opt = epwpl::wpl_ep(data_mat, weight_mat_fit,  
                           rep(v_noise, n_pip), 
                           rep(v_slab, n_pip), 
                           p_incl_grid,
                           damping = .9, k = .99,
                           opt=TRUE, verbose=TRUE,
                           opt_method="Nelder-Mead")
tictoc::toc()

tictoc::tic()
result_vb_opt = epwpl::wpl_vb(data_mat, weight_mat_fit,
                              rep(v_noise, n_pip), 
                              rep(v_slab, n_pip), 
                              p_incl_grid, opt=TRUE)
tictoc::toc()

tictoc::tic()
result_varbvs_opt = epwpl::wpl_varbvs(data_mat, weight_mat_fit,
                                      rep(v_noise, n_pip), 
                                      rep(v_slab, n_pip), 
                                      p_incl_grid, opt=TRUE, verbose=FALSE)
tictoc::tic()



# fix inclusion probability, optimize over variances (the old setting)
result_opt2 = epwpl::wpl_ep(data_mat, weight_mat_fit,
                            v_noise, v_slab, .5,
                            damping = .9, k = .99,
                            opt=TRUE, verbose=TRUE)

# debug individual 2, dim 2, hyperparam 3 (p = .5), singular matrix
i = 2
resp_idx = 2

sqrt_weight = sqrt(weight_mat_fit)
y = data_mat[, resp_idx]
X = data_mat[, -resp_idx]
y_weighted = y * sqrt_weight[i, ]
X_weighted = X * sqrt_weight[i, ]

fit = epwpl::ep_grid_gss(X_weighted, y_weighted, 
                         1, 3, qlogis(0.5),
                         verbose=TRUE, opt = TRUE)

# old EP regression
old_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                               v_noise_grid, v_slab_grid, p_incl_grid,
                               woodbury = FALSE, opt = TRUE)

plot(as.vector(old_result$graphs[[1]]), as.vector(result$graphs[[1]]))

fit = epwpl::ep_grid_gss(X_weighted, y_weighted, 
                         v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                         verbose=FALSE)

# a big grid


# discrete dependent ----

set.seed(1)
n = 100
n_sim = 50
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5 # For Z[i]=-0.1
Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5

Prec1 = Lam1 %*% t(Lam1) + diag(rep(10, p+1))
Prec2 = Lam2 %*% t(Lam2) + diag(rep(10, p+1))

Var1 = solve(Prec1)
Var2 = solve(Prec2)

# covariate matrix
Z = matrix(-.1*(1:n <= n/2)  + .1*(1:n > n/2), nrow = n, ncol = p, byrow = FALSE)

# compute weights
tau = .1  # bandwidth
weight_mat = epwpl::weight_matrix(n, Z, tau)
weight_mat_fit = weight_mat[c(1, n), ]

# true graphs
true_graph_neg = Prec1 != 0
diag(true_graph_neg) = 0
true_graph_pos = Prec2 != 0
diag(true_graph_pos) = 0

X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)
data_mat = rbind(X1, X2)

# initial hyperparameter values
v_noise= 1
p0 = .5
v_slab = 3 

old_result = wpl_ep_regression(data_mat, weight_mat_fit, 
                               v_noise, v_slab, p0,
                               woodbury = FALSE, opt = TRUE)

score_model("EP_opt", epwpl::mean_symmetrize(old_result$graphs[[1]]), 
            true_graph_neg, 1, 1, -.1, p)


# continuous ----
n = 180
p = 4

Prec_cont = function(z) {
  
  STR = 1
  
  pr = matrix(0, p+1, p+1)
  diag(pr) = 2
  
  pr[2,3] = STR
  pr[1,2] = STR*((z>-1) && (z< -.33)) + (STR - STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (0)*((z>0.43) && (z<1))
  pr[1,3] = 0*((z>-1) && (z< -.33)) + (STR*((z+.23)/.56)) * ((z>-0.23) && (z<0.33)) + (STR)*((z>0.43) && (z<1))
  
  pr[2,1] = pr[1,2]
  pr[3,1] = pr[1,3]
  pr[3,2] = pr[2,3]
  
  return(pr)
  # Var = solve(pr)
  # return(Var)
}

Z = c(seq(-0.99, -0.331, (-.331+.99)/59), 
      seq(-0.229,0.329,(.329+.229)/59),
      seq(0.431,.99,(.99-.431)/59))
Z = matrix(Z, n, 1)
X = matrix(0, n, p+1)

p_incl_grid = seq(.05, .8, .05)
v_slab = 3
sigma0 = 1
tau = 0.56

n_sim = 50

# simulate data
for(i in 1:n) {
  X[i, ] = MASS::mvrnorm(1, rep(0, p+1), solve(Var_cont(Z[i])))
}
data_mat = X

# weight matrix
weight_mat = weight_matrix(Z, tau)

tictoc::tic()
ep_vopt_result = 
  epwpl::wpl_ep(data_mat, 
                covariates = Z, tau = tau, weight_mat = weight_mat,
                v_noise_grid = v_noise,
                v_slab_grid = v_slab,
                p_incl_grid = p_incl_grid,
                damping = .9, k = .99,
                opt = TRUE, 
                opt_method = "Nelder-Mead",
                verbose = TRUE)
tictoc::toc()

ep_opt_probs = rbindlist(lapply(1:n, function(i) {
  graph = ep_vopt_result$individuals[[i]]$graph
  incl_prob_model("EP_opt", mean_symmetrize(graph), i, 1, Z[i, ], p) 
}))

ggplot()


