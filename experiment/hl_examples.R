# Small example ----------------------------------------------------------
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
result_epss2 = epwpl::ep_ss2(X, y, v_noise, v_slab, p0,
                             damping=.9, k=.99, opt=TRUE,
                             opt_method = "Nelder-Mead")
# result_epss2 = epwpl::ep_ss2(X, y, v_noise, v_slab, p0,
#                              damping=.9, k=.99, opt=TRUE,
#                              opt_method = "L-BFGS-B",
#                              woodbury = TRUE)
tictoc::toc()


tictoc::tic()
result_epss1 = epwpl::ep_ss1(X, y, v_noise, v_slab, p0, 
                             damping=.9, k=.99, opt=TRUE)
tictoc::toc()

# group spike slab
tictoc::tic()
gss_result = epwpl::GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p0, ncol(X)),
                                      v1 = v_slab, verbose=FALSE, opt=TRUE,
                                      damping = .9, k=.99)
tictoc::toc()

# batch vb
tictoc::tic()
result_vb = epwpl::vb_ss(X, y, v_noise, v_slab, p0)
tictoc::toc()

# varbvs
tictoc::tic()
result_varbvs = varbvs::varbvs(X = X, y = y, Z = NULL, family = "gaussian", 
                               sigma = v_noise, sa = v_slab, logodds = qlogis(p0),
                               update.sigma = TRUE, update.sa = TRUE)
tictoc::toc()

epwpl::normalizelogweights(c(result_varbvs$logw, result_vb$elbo))

## ep with grid prior and importance sampling ----
n_grid = 50
v_noise_grid = rep(v_noise, n_grid)
v_slab_grid = rep(v_slab, n_grid)
p_incl_grid = seq(.01, .9, length.out = n_grid)

tictoc::tic()
bvs_ep2_result = epwpl::epbvs(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt = FALSE, damping = .9, k=1, 
                              method = "ss_ep2")
tictoc::toc()

tictoc::tic()
bvs_ep1_result = epwpl::epbvs(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt = FALSE, damping = .9, k=1, 
                              method = "ss_ep1")
tictoc::toc()

tictoc::tic()
bvs_gss_result = epwpl::epvbs(X,y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                                     opt=FALSE, damping=.9, k=1,
                                    method="gss")
tictoc::toc()

par(mfrow = c(1, 3))
plot(p_incl_grid, bvs_ep2_result$mliks, main="EP 2")
plot(p_incl_grid, bvs_ep1_result$mliks, main="EP 1")
plot(p_incl_grid, bvs_gss_result$mliks, main="GSS")

par(mfrow = c(1, 3))
plot(p_incl_grid, bvs_ep2_result$weights)
plot(p_incl_grid, bvs_ep1_result$weights)
plot(p_incl_grid, bvs_gss_result$weights)

par(mfrow = c(1, 3))
plot(p_incl_grid, bvs_ep2_result$iters, main="EP 2")
plot(p_incl_grid, bvs_ep1_result$iters, main="EP 1")
plot(p_incl_grid, bvs_gss_result$iters, main="GSS")

# compare gss with ss:
pairs(data.frame(
  "EP 2" = log(bvs_ep2_result$pip),
  "EP 1" = log(bvs_ep1_result$pip),
  "GSS"  = log(bvs_gss_result$pip)
))

plot(bvs_ep2_result$pip, gss_grid_result$pip)

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

# my VB
vb_result = epwpl::vb_ss(X, y, v_noise, v_slab, p0, opt=TRUE)
vb_grid_result = vb_grid_ss(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid))

plot(p_incl_grid, vb_grid_result$mliks)

# varbvs
varbvs_result = varbvs::varbvs(X=X, Z=NULL, y=y, family="gaussian", 
                               sigma = v_noise, sa = v_slab, logodds = qlogis(p0))

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

## optimize variance ----
tictoc::tic()
bvs_ep2_result = epwpl::epvbs(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt = TRUE, damping = .9, k=1, 
                              method = "ss_ep2", opt_method = "L-BFGS-B")
tictoc::toc()

tictoc::tic()
bvs_ep1_result = epwpl::epvbs(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt = TRUE, damping = .9, k=1, 
                              method = "ss_ep1")
tictoc::toc()

tictoc::tic()
bvs_gss_result = epwpl::epvbs(X,y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt=TRUE, damping=.9, k=1,
                              method="gss")
tictoc::toc()

par(mfrow = c(1, 3))
plot(p_incl_grid, bvs_ep2_result$mliks, main="EP 2")
plot(p_incl_grid, bvs_ep1_result$mliks, main="EP 1")
plot(p_incl_grid, bvs_gss_result$mliks, main="GSS")




# Bigger example ----------------------------------------------------------

set.seed(1)
d = 100
n = 100
X = cbind(mvtnorm::rmvnorm(n, rep(0, d), 5.*diag(d) + .5*rep(1, d) %*% t(rep(1, d))))
w = c(rep(1, 10), rep(0, d-10))

v_noise = 1/10
y = rnorm(n, X %*% w, sqrt(v_noise))

p0 = .1
v_slab = .1

tictoc::tic()
result = epwpl::ep_ss2(X, y, v_noise, v_slab, .1,
                       opt=TRUE,
                       woodbury = FALSE,
                       opt_method = "UOBYQA")
# result = epwpl::ep_ss2(X, y, v_noise, v_slab, .1, 
#                        opt=TRUE, 
#                        woodbury = FALSE,
#                        opt_method = "Nelder-Mead")
tictoc::toc()

result

tictoc::tic()
gss_result = epwpl::GroupSpikeAndSlab(X, y, 
                                      tau = 1 / v_noise, 
                                      groups = 1:d,
                                      p1 = rep(p0, d),
                                      v1 = v_slab, 
                                      verbose = TRUE, 
                                      opt=TRUE)
tictoc::toc()

tictoc::tic()
result_vb = epwpl::vb_ss(X, y, v_noise, v_slab, p0, opt=TRUE)
tictoc::toc()

## ep with grid prior and importance sampling ----
n_grid = 50
v_noise_grid = rep(v_noise, n_grid)
v_slab_grid = rep(v_slab, n_grid)
p_incl_grid = seq(.1, .9, length.out = n_grid)

tictoc::tic()
bvs_ep2_result = epwpl::epbvs(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt = FALSE, opt_method = "Nelder-Mead", damping = .9, k=.99, 
                              method = "ss_ep2",
                              woodbury = FALSE)
tictoc::toc()

plot(p_incl_grid, bvs_ep2_result$mliks)

tictoc::tic()
bvs_ep1_result = epwpl::epvbs(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt = FALSE, damping = .9, k=1, 
                              method = "ss_ep1")
tictoc::toc()

tictoc::tic()
bvs_gss_result = epwpl::epbvs(X,y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                              opt=FALSE, damping=.9, k=1,
                              method="gss")
tictoc::toc()

plot(p_incl_grid, bvs_gss_result$mliks)

n_grid = 20
v_noise_grid = rep(result$v_noise, n_grid)
v_slab_grid = rep(result$v_slab, n_grid)
p_incl_grid = seq(.1, .5, length.out = n_grid)

result_grid = epwpl::ep_grid_ss(X, y, v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                                opt=FALSE, eps = .9, k= .99)

vb_result_grid = epwpl::vb_grid_ss(X, y, 
                                   v_noise_grid, v_slab_grid, qlogis(p_incl_grid),
                                   opt=TRUE)

varbvs_result_grid = varbvs::varbvs(X=X, Z=NULL, y=y, family="gaussian",
                                    sigma=v_noise_grid, sa=v_slab_grid,
                                    logodds = qlogis(p_incl_grid))

plot(p_incl_grid, bvs_ep2_result$mliks)

plot(p_incl_grid, vb_result_grid$mliks)

plot(p_incl_grid, varbvs_result_grid$logw)

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

