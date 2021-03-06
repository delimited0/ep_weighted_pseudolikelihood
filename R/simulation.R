library(future.apply)
library(doFuture)
library(data.table)

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

score_model = function(method = "method", graph, true_graph, individual, 
                       simulation, covariate, p, w = 1) {
  
  if (is.null(p)) {
    p = nrow(graph) - 1
  }
  
  est_graph = 1 * (graph > 0.5)
  
  data.table(
    method = method,
    sensitivity = sum(est_graph & true_graph) / sum(true_graph),
    specificity = sum(!est_graph & !true_graph) / sum(!true_graph),
    individual = individual, 
    simulation = simulation,
    covariate = covariate,
    p = p,
    w = w
  )
}

#' @export
score_graphs = function(estimates, true_graphs, 
                        threshold = .5, 
                        symmetrizer = mean_symmetrize) {
  
  n_individuals = length(estimates$individuals)
  p = estimates$p
  
  score_table = rbindlist(lapply(1:n_individuals, function(i) {
    
    est = estimates$individuals[[i]]
    tg = true_graphs[[i]]
    est_graph = 1 * (symmetrizer(est$graph) > threshold)
    
    data.table(
      method = estimates$method_name,
      sensitivity = sum(est_graph & tg) / sum(tg),
      specificity = sum(!est_graph & !tg) / sum(!tg),
      individual = i,
      covariate = est$covariates,
      p = p,
      weight1 = 1
    )
  }))
  
  return(score_table) 
}

#' assumptions: estimates1 and estimates2 have same number of individuals, 
#' covariates, p
#' @export
score_combo_graphs = function(estimates1, estimates2, w1, 
                              true_graphs, threshold = .5,
                              symmetrizer = mean_symmetrize) {

  n_individuals = length(estimates1$individuals)
  p = estimates1$p
  
  score_table = rbindlist(lapply(1:n_individuals, function(i) {
    
    est1 = estimates1$individuals[[i]]
    est2 = estimates2$individuals[[i]]
    combo_graph = w1 * est1$graph + (1-w1) * est2$graph
    tg = true_graphs[[i]]
    est_graph = 1 * (symmetrizer(combo_graph) > threshold)
    
    data.table(
      method = paste0(estimates1$method_name, "+", estimates2$method_name),
      sensitivity = sum(est_graph & tg) / sum(tg),
      specificity = sum(!est_graph & !tg) / sum(!tg),
      individual = i,
      covariate = est1$covariates,
      p = p,
      weight1 = w1
    )
  }))
  
  return(score_table)
}

score_likwgt_graphs = function(estimates1, estimates2, 
                               true_graphs, threshold = .5,
                               symmetrizer = mean_symmetrize) {
  n_individuals = length(estimates1$individuals)
  p = estimates1$p
  
  score_table = rbindlist(lapply(1:n_individuals, function(i) {
    
    est1 = estimates1$individuals[[i]]
    est2 = estimates2$individuals[[i]]
    
    lliks = data.table(
      est1 = rowMeans(est1$llik),
      est2 = rowMeans(est2$llik)
    )
    
    
    
    combo_graph = w1 * est1$graph + (1-w1) * est2$graph
    tg = true_graphs[[i]]
    est_graph = 1 * (symmetrizer(combo_graph) > threshold)
    
    data.table(
      method = paste0(estimates1$method_name, "+", estimates2$method_name),
      sensitivity = sum(est_graph & tg) / sum(tg),
      specificity = sum(!est_graph & !tg) / sum(!tg),
      individual = i,
      covariate = est1$covariates,
      p = p,
      weight1 = w1
    )
  }))
  
  return(score_table)
}

incl_prob_model = function(method = "method", graph, individual, simulation, 
                           covariate, p, w1 = 1) {
  
  data.table(
    method = method,
    incl_prob_12 = graph[1, 2],
    incl_prob_13 = graph[1, 3],
    individual = individual,
    simulation = simulation,
    covariate = covariate,
    p = p,
    weight1 = w1
  )
}

#' @export
incl_prob_123 = function(estimates, symmetrizer = mean_symmetrize) {
  
  n_individuals = length(estimates$individuals)
  p = estimates$p
  
  ip_table = rbindlist(lapply(1:n_individuals, function(i) {
    est = estimates$individuals[[i]]
    est_graph = symmetrizer(est$graph)
    # est_graph = est$graph
    
    data.table(
      method = estimates$method_name,
      incl_prob_12 = est_graph[1, 2],
      incl_prob_13 = est_graph[1, 3],
      individual = i,
      covariate = est$covariates,
      p = p,
      weight1 = 1
    )
  }))
  
  return(ip_table)
}

#' @export
incl_prob_123_combo = function(estimates1, estimates2, w1,
                               symmetrizer = mean_symmetrize) {
  
  n_individuals = length(estimates1$individuals)
  p = estimates1$p
  
  score_table = rbindlist(lapply(1:n_individuals, function(i) {
    
    est1 = estimates1$individuals[[i]]
    est2 = estimates2$individuals[[i]]
    combo_graph = w1 * est1$graph + (1-w1) * est2$graph
    est_graph = symmetrizer(combo_graph)
    
    data.table(
      method = paste0(w1, "_", estimates1$method_name, "+", 
                      1-w1, "_", estimates2$method_name),
      incl_prob_12 = est_graph[1, 2],
      incl_prob_13 = est_graph[1, 3],
      individual = i,
      covariate = est1$covariates,
      p = p,
      weight1 = w1
    )
  }))
  
  return(score_table)
}

alpha_123 = function(estimates, symmetrizer = mean_symmetrize) {
  n_individuals = length(estimates$individuals)
  p = estimates$p
  
  ip_table = rbindlist(lapply(1:n_individuals, function(i) {
    
    est = estimates$individuals[[i]]
    
    data.table(
      method = estimates$method_name,
      alpha_12 = est$fits[[2]]$alpha[1, ],
      alpha_13 = est$fits[[3]]$alpha[1, ],
      w_12 = est$fits[[2]]$w,
      w_13 = est$fits[[3]]$w,
      logodds_12 = est$fits[[2]]$logodds,
      logodds_13 = est$fits[[3]]$logodds,
      individual = i,
      covariate = est$covariates,
      p = p,
      weight1 = 1
    )
  }))
}

#' @export
weight_matrix = function(covariates, tau) {
  
  n = nrow(covariates)
  p = ncol(covariates)
  
  weight_mat = matrix(1, n, n)
  for(i in 1:n){
    for(j in 1:n){
      weight_mat[i, j] = dnorm(norm(covariates[i, ] - covariates[j, ], "2"), 0, tau)
    }
  }
  
  for(i in 1:n){
    weight_mat[, i] = n * (weight_mat[, i] / sum(weight_mat[, i])) #Scaling the weights so that they add up to n
  }
  
  return(weight_mat)
}

vsvb_to_graphs = function(vsvb_list) {
  n = nrow(vsvb_list[[1]])
  p = length(vsvb_list)
  graphs = replicate(n, diag(p), simplify = FALSE)
  
  for (dimen in 1:p) {
    for (i in 1:n) {
      prob_row = rep(0, p)
      prob_row[-dimen] = vsvb_list[[dimen]][i, ]
        
      graphs[[i]][dimen, ] = prob_row
    }  
  }
  
  return(graphs)
}
