# My own example ----------------------------------------------------------
set.seed(1)

d = 10
n = 100
X = matrix(rnorm(n * d), nrow = n, ncol=d)

v_noise = .1
v_slab = 1
p0 = .2
# z = rbinom(d, 1, p0)
# w = rep(0, d)
# w[z == 1] = rnorm(sum(z == 1), 0, sqrt(v_slab))
# w = c(1, rep(0, d-1))
w = c(1, 1, rep(0, d-2))

y = rnorm(n, X %*% w, sqrt(v_noise))

# first dim plot
par(mfrow = c(2, 2))
plot(X[, 1], y)
plot(X[, 2], y)
plot(X[, 3], y)
plot(X[, 4], y)

tictoc::tic()
result = epwpl::ep_wlr(X, y, v_noise, v_slab, p0, max_iter = 1000)
tictoc::toc()
result$m

# group spike slab
tictoc::tic()
gss_result = epwpl::GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p0, ncol(X)),
                                      v1 = v_slab, verbose=FALSE, opt=TRUE)
tictoc::toc()

# ep with grid prior and importance sampling
n_grid = 50
v_noise_grid = rep(v_noise, n_grid)
v_slab_grid = rep(v_slab, n_grid)
p_incl_grid = seq(.01, .9, length.out = n_grid)
# p_incl_grid = seq(.1, .12, length.out = n_grid)
# p_incl_grid = c(.109, .110, .111, .112)
#
is_grid_result = epwpl::ep_grid_lr(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                           opt = FALSE, eps = 1, k=.99)

par(mfrow = c(1, 1))
plot(p_incl_grid, is_grid_result$mliks)
plot(p_incl_grid, is_grid_result$weights)

# ep gss with prior grid and inportance sampling
n_grid = 10
v_noise_grid = rep(v_noise, n_grid)
v_slab_grid = rep(v_slab, n_grid)
p_incl_grid = seq(.01, .9, length.out = n_grid)

tictoc::tic()
gss_grid_result = epwpl::ep_grid_gss(X,y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                                opt=FALSE)
tictoc::toc()

plot(p_incl_grid, gss_grid_result$mliks)
plot(p_incl_grid, gss_grid_result$weights)

# compare gss with ss:
plot(is_grid_result$pip, gss_grid_result$pip)

# the posteriors are the same
plot(is_result$alpha[, 1], gss_result$alpha[, 1])
plot(is_result$alpha[, 60], gss_result$alpha[, 60])
plot(is_result$alpha[, n_grid], gss_result$alpha[, n_grid])

#
plot(p_incl_grid, gss_result$weights, col="blue")
points(p_incl_grid, is_result$weights)

# compare to carbonetto stephens
varbvs_result = varbvs::varbvs(X = X, y = y, Z = NULL) 
                               # sigma = v_noise_grid, 
                               # sa = v_slab_grid,
                               # logodds = qlogis(p_incl_grid))

plot(p_incl_grid, varbvs_result$logw)

plot(is_result$sigma, varbvs_result$sigma)
plot(is_result$sa, varbvs_result$sa)

plot(is_result$alpha)

# ep, also optimize p0
result_p0 = epwpl::ep_wlr_nmp(X, y, sigma0, p0, v_slab)

microbenchmark::microbenchmark(
  result = epwpl::ep_wlr(X, y, sigma0, p0, v_slab, max_iter = 2000)
)

# compare to mcmc
prior = BoomSpikeSlab::SpikeSlabPrior(X, y, expected.model.size = 1)
mcmc_result = BoomSpikeSlab::lm.spike(y ~ X - 1, niter = 5000)

# compare to home rolled vsvb
dp_result = vb_wlr(X, y, result$sigma0, p0, result$v_slab)
dp_result2 = vb_wlr(X, y, result$sigma0, .9, result$v_slab)



# Bigger example ----------------------------------------------------------

d = 100
n = 100
X = cbind(mvtnorm::rmvnorm(n, rep(0, d), 5.*diag(d) + .5*rep(1, d) %*% t(rep(1, d))))
w = c(rep(1, 10), rep(0, d-10))

v_noise = 1/10
y = rnorm(n, X %*% w, sqrt(v_noise))

p0 = .1
v_slab = .1

tictoc::tic()
result = epwpl::ep_wlr(X, y, v_noise, v_slab, p0, 
                       max_iter = 1000,
                       opt=TRUE,
                       woodbury = FALSE)
tictoc::toc()

tictoc::tic()
gss_result = epwpl::GroupSpikeAndSlab(X, y, 
                                      tau = 1 / v_noise, 
                                      groups = 1:d,
                                      p1 = rep(p0, d),
                                      v1 = v_slab, 
                                      verbose = TRUE, 
                                      opt=TRUE)
tictoc::toc()

n_grid = 20
v_noise_grid = rep(result$v_noise, n_grid)
v_slab_grid = rep(result$v_slab, n_grid)
p_incl_grid = seq(.1, .5, length.out = n_grid)

result_grid = epwpl::ep_grid_ss(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                                opt=FALSE, eps = .9, k= .99)

plot(p_incl_grid, result_grid$mliks)

# Example 5.1 -------------------------------------------------------------

d = 2
n = 2
X = mvtnorm::rmvnorm(n, rep(0, d), 5.*diag(d) + .5*rep(1, d) %*% t(rep(1, d)))

p0 = .5
v_slab = 1
z = rbinom(d, 1, p0)
w = rep(0, d)
w[z == 1] = rnorm(sum(z == 1), 0, sqrt(v_slab))

sigma_noise = sqrt(1/10)
y = rnorm(n, X %*% w, sigma_noise)

weights = rep(1, n)
sigma0 = sqrt(.1)
result = ep_wlr(X, y, weights, sigma0, p0, v_slab, max_iter = 2000)

foo = matrix(1:9, nrow = 3)

foo * c(1, 0, -1)
foo %*% diag(c(1, 0, -1))

# compare with varbvs on its doc example ----------------------------------
maf <- 0.05 + 0.45*runif(200)
X   <- (runif(400*200) < maf) + (runif(400*200) < maf)
X   <- matrix(as.double(X),400,200,byrow = TRUE)
# Z   <- randn(400,3)

# u    <- c(-1,2,1)
beta <- c(rnorm(20),rep(0,180))
beta <- 1/sd(c(X %*% beta)) * beta
# y <- c(-2 + Z %*% u + X %*% beta + rnorm(400))
y <- c(-2 + X %*% beta + rnorm(400))

logodds_grid = seq(-3,-1,0.1)
sigma_grid = rep(1, length(logodds_grid))
sa_grid = rep(1, length(logodds_grid))

ep_fit = epwpl::ep_grid_lr(X, y, sigma_grid, sa_grid, logodds_grid, opt=FALSE)

fit <- varbvs::varbvs(X, Z = NULL, y, logodds = logodds_grid,
                      update.sigma=FALSE, update.sa=FALSE)

plot(plogis(logodds_grid), ep_fit$weights)

plot(ep_fit$pip, fit$pip)

# compare with varbvs on a small example -----------------------------------
maf <- 0.05 + 0.45*runif(20)
X   <- (runif(40*20) < maf) + (runif(400*200) < maf)
X   <- matrix(as.double(X),400,200,byrow = TRUE)
# Z   <- randn(400,3)

# u    <- c(-1,2,1)
beta <- c(rnorm(20),rep(0,180))
beta <- 1/sd(c(X %*% beta)) * beta
# y <- c(-2 + Z %*% u + X %*% beta + rnorm(400))
y <- c(-2 + X %*% beta + rnorm(400))

