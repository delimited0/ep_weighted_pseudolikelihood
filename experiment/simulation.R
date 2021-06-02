library(future.apply)
library(doFuture)
library(data.table)

wpl_regression = function(data_mat, weight_mat, sigma0, p0, v_slab, 
                          n_threads = 1, blas_threads = 1, woodbury = FALSE) {
  registerDoFuture()
  plan(multisession, workers = n_threads)
  RhpcBLASctl::blas_set_num_threads(blas_threads)
  RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  n = nrow(weight_mat)
  
  sqrt_weight = sqrt(weight_mat)

  progressr::with_progress({
    prog = progressr::progressor(along = 1:n)

    graphs = future_lapply(1:n, function(i) {

      prog(sprintf("Individual =%g", i))

      incl_prob = sapply(1:p, function(resp_idx) {

        y = data_mat[, resp_idx]
        X = data_mat[, -resp_idx]
        y_weighted = y * sqrt_weight[i, ]
        X_weighted = X * sqrt_weight[i, ]

        fit = epwpl::ep_wlr(X_weighted, y_weighted, sigma0, p0, v_slab,
                            woodbury = woodbury)

        prob_row = matrix(0, nrow = 1, ncol = p)
        prob_row[, -resp_idx] = t(plogis(fit$p))
        return(prob_row)
      })

      return(incl_prob)
    })
  })

  return(graphs)
}

mean_symmetrize = function(mat) {
  for(i in 1:(p+1)) {
    for(j in i:(p+1)) {
      mat[i, j] = mean(c(mat[i, j], mat[j, i]))
      mat[j, i] = mat[i, j]
    }
  }
  return(mat)
}

max_symmetrize = function(mat) {
  for(i in 1:(p+1)) {
    for(j in i:(p+1)) {
      mat[i, j] = max(mat[i, j], mat[j, i])
      mat[j, i] = mat[i, j]
    }
  }
  return(mat)
}

score_model = function(graph, true_graph, individual, simulation, covariate, p) {
  
  if (is.null(p)) {
    p = nrow(graph) - 1
  }
  
  est_graph = 1 * (graph > 0.5)
  
  data.table(
    sensitivity = sum(est_graph & true_graph) / sum(true_graph),
    specificity = sum(!est_graph & !true_graph) / sum(!true_graph),
    individual = individual, 
    simulation = simulation,
    covariate = covariate,
    p = p
  )
}

weight_matrix = function(n, cov_mat) {
  p = ncol(cov_mat)
  
  weight_mat = matrix(1, n, n)
  for(i in 1:n){
    for(j in 1:n){
      weight_mat[i, j] = dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
    }
  }
  
  for(i in 1:n){
    weight_mat[, i] = n * (weight_mat[, i] / sum(weight_mat[, i])) #Scaling the weights so that they add up to n
  }
  
  return(weight_mat)
}