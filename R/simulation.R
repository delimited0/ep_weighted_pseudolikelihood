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
    est_graph = symmetrizer(1 * (est$graph > threshold))
    
    data.frame(
      method = estimates$method_name,
      sensitivity = sum(est_graph & tg) / sum(tg),
      specificity = sum(!est_graph & !true_graph) / sum(!true_graph),
      individual = i,
      covariate = est$covariates,
      p = p
    )
  }))
  
  return(score_table) 
}

incl_prob_model = function(method = "method", graph, individual, simulation, 
                           covariate, p) {
  
  data.table(
    method = method,
    incl_prob_12 = graph[1, 2],
    incl_prob_13 = graph[1, 3],
    individual = individual,
    simulation = simulation,
    covariate = covariate,
    p = p
  )
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
