library(future.apply)
library(doFuture)
library(data.table)

wpl_ep_regression = function(data_mat, weight_mat, sigma0, p0, v_slab, 
                          woodbury = FALSE, opt = TRUE) {
  # registerDoFuture()
  # plan(multisession, workers = n_threads)
  # RhpcBLASctl::blas_set_num_threads(blas_threads)
  # RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  p = ncol(data_mat)
  n = nrow(weight_mat)
  
  sqrt_weight = sqrt(weight_mat)
  
  estimated_sigma_noise = rep(NA, p)
  estimated_v_slab = rep(NA, p)
  graphs = replicate(n, matrix(0, p, p), simplify = FALSE)

  progressr::with_progress({
    prog = progressr::progressor(along = 1:n)
    
    for (i in 1:n) {
      prog(sprintf("Individual =%g", i))
      
      for (resp_idx in 1:p) {
        
        y = data_mat[, resp_idx]
        X = data_mat[, -resp_idx]
        y_weighted = y * sqrt_weight[i, ]
        X_weighted = X * sqrt_weight[i, ]
        
        fit = epwpl::ep_wlr(X_weighted, y_weighted, sigma0, p0, v_slab,
                            woodbury = woodbury, opt = opt)
        
        graphs[[i]][resp_idx, -resp_idx] = t(plogis(fit$p))
        estimated_sigma_noise[resp_idx] = fit$sigma0
        estimated_v_slab[resp_idx] = fit$v_slab
      }
    }
  })

  return(list(graphs = graphs,
              sigma_noise = estimated_sigma_noise,
              v_slab = estimated_v_slab))
}

wpl_vb_regression = function(data_mat, weight_mat, sigma0, p0, v_slab) {
  
  p = ncol(data_mat)
  n = nrow(weight_mat)
  
  if (length(sigma0) == 1) {
    sigma0 = rep(sigma0, p)
  }
  if (length(v_slab) == 1) {
    v_slab = rep(v_slab, p)
  }
  
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
        
        browser()
        fit = epwpl::vb_wlr(X_weighted, y_weighted, 
                            sigma0[resp_idx], p0, v_slab[resp_idx])
        
        prob_row = matrix(0, nrow = 1, ncol = p)
        prob_row[, -resp_idx] = t(plogis(fit$logit_alpha))
        return(prob_row)
      })
      
      return(incl_prob)
    })
  })
  
  return(graphs)
}

wpl_vsvb_regression = function(data_mat, weight_mat, sigma0, p0, v_slab,
                               tune = FALSE, blas_threads = 1) {
  
  RhpcBLASctl::blas_set_num_threads(blas_threads)
  RhpcBLASctl::omp_set_num_threads(blas_threads)
  
  n = nrow(data_mat)
  p = ncol(data_mat)-1
  
  mylist = rep(list(matrix(0, n, p)), p+1)
  
  for (resp_index in 1:(p+1)) {
    
    y = data_mat[, resp_index] #Set variable number `resp_index` as the response
    X_mat = data_mat[, -resp_index]
    X_vec = matrix(0, n*p, 1)
    X = matrix(rep(0, n^2*p), nrow = n, ncol = n*p)
    
    for(i in 1:n) {
      for(j in 1:p) {
        k = p*(i-1) + 1 
        X[i, k+j-1] = X_mat[i, j]
        X_vec[k+j-1] = X[i, k+j-1]
      }
    }
    
    Big_diag_mat <- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)
    for(i in 1:n){
      k=p*(i-1)
      for(j in 1:p){
        Big_diag_mat[i,k+j]=1
      }
    }
    
    sigmasq = sigma0^2
    
    XtX = t(X) %*% X
    
    DXtX = diag(XtX)
    DXtX_rep = rep(DXtX, p) 
    DXtX_mat = matrix(DXtX_rep, n*p, p, byrow=FALSE)
    
    alpha = rep(.2, n*p)
    sigmabeta_sq = v_slab
    mu = rep(0, p) # Variational parameter
    true_pi = p0
    
    y_long_vec = as.vector(t(y %*% matrix(1, 1, p)))
    Xty=t(X) %*% as.vector(y)
    mu_mat = matrix(NA, n, p)
    
    D_long = matrix(0, n*p, n)
    for( i in 1:n){
      D_long[, i] = matrix(t(weight_mat[, i] %*% matrix(1, 1, p)), n*p, 1)
    }
    
    S_sq = matrix(sigmasq*(DXtX + 1/sigmabeta_sq)^(-1), n, p) #Initialization
    
    iter = 1
    
    ind_vec = seq(0,(n-1)*p,by=p)
    Ind_mat = matrix(0,n,p)
    for(j in 1:p){
      Ind_mat[, j] = ind_vec + j
    }
    Big_ind = matrix(0, n*p, p)
    Big_ind_1 = matrix(0, n*p, p)
    for(j in 1:p){
      Big_ind[Ind_mat[, j], j] = 1
      Big_ind_1[Ind_mat[, j], j] = 0
    }
    
    DXtX_Big_ind = DXtX_mat * Big_ind
    
    candL = seq(0.1, 0.9, .2)# Different values of hyperparameter true_pi
    #candL=0.5
    like = rep(0, length(candL))
    elb = like
    
    est_pi = rep(0, n)
    est_q = est_pi
    
    ####################tuning hyperparameters##################################
    if (tune) {
      idmod = varbvs::varbvs(X_mat, y, Z=Z[, 1], verbose = FALSE)#Setting hyperparameter value as in Carbonetto Stephens model
      inprob = idmod$pip
      rest_index_set = setdiff(c(1:(p+1)), resp_index)
      
      sigmasq = mean(idmod$sigma)
      pi_est = mean(1 / (1 + exp(-idmod$logodds)))
      sigmavec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)
      elb1 = matrix(0, length(sigmavec), 1)
      
      for (j in 1:length(sigmavec)) {
        res = epwpl::cov_vsvb(y, X, X_mat, mu, mu_mat, alpha, DXtX_Big_ind, 
                              weight_mat, D_long, sigmasq, sigmabeta_sq,
                              y_long_vec, X_vec, true_pi)
        elb1[j] = res$var.elbo
        
      }
      sigmabeta_sq = sigmavec[which.max(elb1)] #Choosing hyperparameter based on ELBO maximization
      true_pi = pi_est
    }
    
    result = epwpl::cov_vsvb(y, X, X_mat, mu, mu_mat, alpha, DXtX_Big_ind, 
                             weight_mat, D_long, sigmasq, sigmabeta_sq,
                             y_long_vec, X_vec, true_pi)
    
    incl_prob = result$var.alpha
    mu0_val = result$var.mu0_lambda
    
    heat_alpha = matrix(incl_prob, n, p, byrow=TRUE)
    
    mylist[[resp_index]] = heat_alpha
  }
  
  result = list(
    graphs = vsvb_to_graphs(mylist),
    sigma0sq = sigmasq,
    v_slab = sigmabeta_sq
  )
  
  return(result)
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

score_model = function(method = "method", graph, true_graph, individual, 
                       simulation, covariate, p) {
  
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
