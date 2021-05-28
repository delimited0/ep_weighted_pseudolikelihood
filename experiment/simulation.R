library(future.apply)
library(doFuture)

wpl_regression = function(data_mat, weight_mat, sigma0, p0, v_slab, 
                          n_threads = 1, blas_threads = 1, woodbury = FALSE) {
  # registerDoFuture()
  # plan(multisession, workers = n_threads)
  RhpcBLASctl::blas_set_num_threads(blas_threads)
  RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  n = nrow(data_mat)
  
  sqrt_weight = sqrt(weight_mat)

  progressr::with_progress({
    prog = progressr::progressor(along = 1:n)

    graphs = lapply(1:n, function(i) {

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

  # progressr::with_progress({
  #   
  #   prog = progressr::progressor(along = 1:n*p)
  #   
  #   regressions =
  #     foreach(i = 1:n) %:%
  #       foreach(resp_idx = 1:p, .combine = 'rbind') %dopar% {
  #       
  #         prog(sprintf("Ind=%g, Cov=%g", i, resp_idx))
  #         
  #         y = data_mat[, resp_idx]
  #         X = data_mat[, -resp_idx]
  #         y_weighted = y * sqrt_weight[i, ]
  #         X_weighted = X * sqrt_weight[i, ]
  #         
  #         fit = epwpl::ep_wlr(X_weighted, y_weighted, sigma0, p0, v_slab)
  #         
  #         prob_row = matrix(0, nrow = 1, ncol = p)
  #         prob_row[, -resp_idx] = t(plogis(fit$p))
  #         return(prob_row)
  #       }
  # })
  # 
  return(graphs)
}
