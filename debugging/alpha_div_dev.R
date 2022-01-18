# tiny example ----
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

alpha = .5
result_alpha = epwpl::alpha_ss(X, y, alpha, v_noise, v_slab, p0, 
                               n_div = 10, lr = .01, min_logval = -500)

result_l0 = epwpl::l0_ss(X, y, v_noise, v_slab, p0)


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

alpha = .5
result_alpha = epwpl::alpha_ss(X, y, alpha, v_noise, v_slab, p0, 
                               n_div = 10, lr = .01, min_logval = -50)
