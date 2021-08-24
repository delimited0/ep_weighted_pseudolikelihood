


# My own example ----------------------------------------------------------

d = 10
n = 100
X = cbind(1, mvtnorm::rmvnorm(n, rep(0, d), 5.*diag(d) + .5*rep(1, d) %*% t(rep(1, d))))

v_slab = .1
# z = rbinom(d, 1, p0)
# z = c(1, rep(0, d))
# w = rep(0, d)
# w[z == 1] = rnorm(sum(z == 1), 0, sqrt(v_slab))
w = c(1, 1, rep(0, d-1))

sigma_noise = sqrt(1/10)
# sigma_noise = .1
y = rnorm(n, X %*% w, sigma_noise)

p0 = .1
weights = rep(1, n)
sigma0 = sqrt(1)

tictoc::tic()
result = epwpl::ep_wlr(X, y, sigma0, p0, v_slab, max_iter = 2000)
tictoc::toc()
result$m

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

# compare to carbonetto stephens
varbvs_result = varbvs::varbvs(X = X, y = y, Z = NULL)

# Bigger example ----------------------------------------------------------

d = 1000
n = 100
X = cbind(1, mvtnorm::rmvnorm(n, rep(0, d), 5.*diag(d) + .5*rep(1, d) %*% t(rep(1, d))))
w = c(1, rep(1, 10), rep(0, d-10))

sigma_noise = sqrt(1/10)
y = rnorm(n, X %*% w, sigma_noise)

p0 = .1
weights = rep(1, n)
sigma0 = sqrt(1)
v_slab = .1

tictoc::tic()
result = epwpl::ep_wlr(X, y, sigma0, p0, v_slab, max_iter = 2000)
tictoc::toc()
result$m

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
