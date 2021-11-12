s <- Rfast::cova( as.matrix(iris[, 1:4]) )

microbenchmark::microbenchmark(
  # Rfast::spdinv(s),
  solve(s),
  chol2inv(chol(s))
  # check = 'equal'
)

determinant(s)



d = 500
d = 3
A = matrix(1:d^2, nrow = d)

microbenchmark::microbenchmark(
  A %*% diag(1:d) %*% t(A),
  t(t(A) * 1:d) %*% t(A),
  sweep(A, MARGIN=2, 1:d, `*`) %*% t(A),
  # A * outer(rep.int(1L, nrow(A)), 1:d) %*% t(A),
  check = 'equal'
)

microbenchmark::microbenchmark(
  diag(1:d) %*% t(A),
  t(A) * 1:d,
  check = 'equal'
)


d = 1000
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))

microbenchmark::microbenchmark(
  determinant(Sigma)$modulus,
  2*sum(log(diag(chol(Sigma)))),
  sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
)
