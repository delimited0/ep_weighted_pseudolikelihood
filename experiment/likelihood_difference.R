
# 2d example --------------------------------------------------------------
set.seed(1)
library(data.table)
library(ggplot2)

d = 2
n = 1
X = matrix(rnorm(n * d), nrow = n, ncol=d)

v_noise = .1
v_slab = 1
p0 = .2

w = c(1, 0)
y = rnorm(n, X %*% w, sqrt(v_noise))

result_epss2 = epwpl::ep_ss2(X, y, v_noise, v_slab, p0, 
                             damping=.9, k=.99, opt=FALSE)
result_varbvs =  varbvs::varbvs(X = X, y = y, Z = NULL, family = "gaussian", 
                                sigma = v_noise, sa = v_slab, logodds = log10(p0 / (1 - p0)),
                                update.sigma = FALSE, update.sa = FALSE)

result_epss2$llik
result_varbvs$logw

# Small example ----------------------------------------------------------
set.seed(1)
library(data.table)
library(ggplot2)

d = 10
n = 100
X = matrix(rnorm(n * d), nrow = n, ncol=d)

v_noise = 1
v_slab = 1
p0 = .2
w = c(1, 1, rep(0, d-2))

y = rnorm(n, X %*% w, sqrt(v_noise))


## no optimization ----

result_epss2 = epwpl::ep_ss2(X, y, v_noise, v_slab, p0, 
                             damping=1, k=1, opt=FALSE, woodbury = TRUE, max_iter = 1)
result_epss1 = epwpl::ep_ss1(X, y, v_noise, v_slab, p0, 
                             damping=.9, k=.99, opt=FALSE)
result_vb = epwpl::vb_ss(X, y, v_noise, v_slab, p0, opt = FALSE)
result_varbvs =  varbvs::varbvs(X = X, y = y, Z = NULL, family = "gaussian", 
                                sigma = v_noise, sa = v_slab, logodds = log10(p0 / (1 - p0)),
                                update.sigma = FALSE, update.sa = FALSE, verbose=FALSE)
result_gss = epwpl::GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p0, ncol(X)),
                                      v1 = v_slab, verbose=FALSE, opt=FALSE,
                                      damping = .9, k=.99)

all_result = rbind(
  data.table(
    "EP" = result_epss2$m[, 1],
    "GSS" = result_gss$meanMarginals[, 1],
    "VB" = result_vb$mu[, 1] * result_vb$alpha[,1],
    "varbvs" = result_varbvs$mu[, 1] * result_varbvs$alpha[,1],
    "Estimate" = "Mean",
    "dim" = 1:d
  ), 
  data.table(
    "EP" = result_epss2$p[,1],
    "GSS" = plogis(result_gss$posteriorApproximation$p),
    "VB" = result_vb$alpha[,1],
    "varbvs" = result_varbvs$alpha[, 1],
    "Estimate" = "p_incl",
    "dim" = 1:d
  ),
  data.table(
    "EP" = result_epss2$v,
    "GSS" = result_gss$varMarginals,
    "VB" = result_vb$v * result_vb$alpha[,1] + 
      result_vb$mu[,1]^2 * (1 - result_vb$alpha[,1]) * result_vb$alpha[,1],
    "varbvs" = result_varbvs$s[, 1] * result_varbvs$alpha[, 1] +
      result_varbvs$mu[,1]^2 * (1 - result_varbvs$alpha[,1]) * result_varbvs$alpha[,1],
    "Estimate" = "Variance",
    "dim" = 1:d
  ),
  data.table(
    "EP" = result_epss2$v_noise,
    "GSS" = result_gss$v_noise,
    "VB" = result_vb$v_noise,
    "varbvs" = result_varbvs$sigma,
    "Estimate" = "v_noise",
    "dim" = 1:d
  ),
  data.table(
    "EP" = result_epss2$v_slab,
    "GSS" = result_gss$v_slab,
    "VB" = result_vb$v_slab,
    "varbvs" = result_varbvs$sa,
    "Estimate" = "v_slab",
    "dim" = 1:d
  )
)

long_result = melt(all_result, id.vars = c("Estimate", "dim"), variable.name = "Method")

ggplot(long_result, aes(x = dim, y = value, color = Method)) +
  geom_point() +
  facet_wrap(vars(Estimate), scales = "free")

result_epss2$llik
result_vb$elbo
result_varbvs$logw
result_gss$evidence
result_epss1$llik

plot(result_epss2$v, result_gss$varMarginals)


## optimization ----

result_epss2 = epwpl::ep_ss2(X, y, v_noise, v_slab, p0, 
                             damping=.9, k=.99, opt=TRUE, 
                             opt_method = "Nelder-Mead")
result_vb = epwpl::vb_ss(X, y, v_noise, v_slab, p0)
result_varbvs = varbvs::varbvs(X = X, y = y, Z = NULL, family = "gaussian", 
                               sigma = v_noise, sa = v_slab, logodds = log10(p0 / (1 - p0)),
                               update.sigma = TRUE, update.sa = TRUE)
prior = BoomSpikeSlab::SpikeSlabPrior(X, y, expected.model.size = 1)
result_mcmc = BoomSpikeSlab::lm.spike(y ~ X - 1, niter = 5000)

# estimates
all_result = rbind(
  data.table(
    "EP" = result_epss2$m[,1],
    "VB" = result_vb$mu[,1],
    "varbvs" = result_varbvs$mu[, 1],
    "Estimate" = "Mean",
    "dim" = 1:d
  ), 
  data.table(
    "EP" = result_epss2$p[,1],
    "VB" = result_vb$alpha[,1],
    "varbvs" = result_varbvs$alpha[, 1],
    "Estimate" = "p_incl",
    "dim" = 1:d
  ),
  data.table(
    "EP" = result_epss2$v,
    "VB" = result_vb$v,
    "varbvs" = result_varbvs$s[, 1],
    "Estimate" = "Variance",
    "dim" = 1:d
  ),
  data.table(
    "EP" = result_epss2$v_noise,
    "VB" = result_vb$v_noise,
    "varbvs" = result_varbvs$sigma,
    "Estimate" = "v_noise",
    "dim" = 1:d
  ),
  data.table(
    "EP" = result_epss2$v_slab,
    "VB" = result_vb$v_slab,
    "varbvs" = result_varbvs$sa,
    "Estimate" = "v_slab",
    "dim" = 1:d
  )
)

library(ggplot2)

long_result = melt(all_result, id.vars = c("Estimate", "dim"), variable.name = "Method")

ggplot(long_result, aes(x = dim, y = value, color = Method)) +
  geom_point() +
  facet_wrap(vars(Estimate), scales = "free")

result_epss2$llik
result_vb$elbo
result_varbvs$logw


# HL toy example ----------------------------------------------------------

mu = c(0, 0)
Sigma = matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2)
v_slab = 1
p0 = .5
v_noise = 1 / 10

# n_sim = 100000
n_train = 100

X = MASS::mvrnorm(n_train, mu, Sigma)
z = rbinom(n_train, 1, p0)
# w = z * rnorm(2, 0, sqrt(v_slab))
w = c(0, 1)
y = rnorm(n_train, X %*% w, sqrt(v_noise))

result_epss2 = epwpl::ep_ss2(X, y, v_noise, v_slab, p0, 
                             damping=.9, k=.99, opt=FALSE)
result_varbvs = varbvs::varbvs(X = X, y = y, Z = NULL, family = "gaussian", 
                               sigma = v_noise, sa = v_slab, logodds = log10(p0 / (1 - p0)),
                               update.sigma = FALSE, update.sa = FALSE, verbose=FALSE)
# result_gss = epwpl::GroupSpikeAndSlab(X, y, tau=1/v_noise, p1 = rep(p0, ncol(X)),
#                                       v1 = v_slab, verbose=FALSE, opt=FALSE,
#                                       damping = .9, k=.99)

prior = BoomSpikeSlab::SpikeSlabPrior(X, y, expected.model.size = 1)
mcmc_result = BoomSpikeSlab::lm.spike(y ~ X - 1, niter = 10000)

library(ellipse)
library(ggplot2)

ep_fit = data.table(ellipse(x = result_epss2$V, centre = result_epss2$m))
vb_fit = data.table(ellipse(x = diag(result_varbvs$s[, 1] * result_varbvs$alpha[, 1] +
                                       result_varbvs$mu[,1]^2 * 
                                       (1 - result_varbvs$alpha[,1]) * result_varbvs$alpha[,1]),
                            centre = result_varbvs$mu * result_varbvs$alpha[, 1]))
# gss_fit = data.table(ellipse(x = result_gss$V, centre = result_gss$m))
ggplot(ep_fit, aes(x = x, y = y)) +
  geom_path(color = "blue") +
  geom_path(data = vb_fit, aes(x = x, y = y), color = "red", lty =2) +
  geom_point(data = data.table(mcmc_result$beta[5000:10000,]), aes(x = X1, y = X2))
  


  





