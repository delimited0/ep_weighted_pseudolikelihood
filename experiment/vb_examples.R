library(ggplot2)
library(data.table)

# Discrete covariate, independent ------------------------------------------------------
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3))*5#For Z[i]=-0.1
Lam2 = Lam1 #Same lambda for both covariate levels, corresponds to covariate independent levels

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)

Z = matrix(-1, n, p) #Initializing the covariate matrix

beta = matrix(0, n, p) # Ground truth of the dependence structure
resp_index = 1;# The index we consider as response
mylist = rep(list(beta), p+1) #The variable specific inclusion probability matrix:ith row corresponds to the dependence structure for the i th subject, j th matrix corresponds to 
# the j th variable as response and the remaining as predictors.
data_mat = rbind(X1, X2)


for (resp_index in 1:(p+1)) {
  
  for (i in 1:n) {
    beta[i, ] = (t(Lam1[-resp_index]) > 0) * (i <= n/2) + 
      (t(Lam2[-resp_index]) > 0) * (i > n/2) #Ground truth
    
    for(j in 1:p){
      Z[i,j] = -.1*(i <= n/2)  + .1*(i > n/2)
    }
  }
  
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
  
  ELBO_LBit=rep(0,10000)
  Big_diag_mat <- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)
  for(i in 1:n){
    k=p*(i-1)
    for(j in 1:p){
      Big_diag_mat[i,k+j]=1
    }
  }
  q = matrix(2, n, 1)
  sigmasq = 1 #Initialization of the hyperparameter value
  E <- rnorm(n, 0, sigmasq)
  
  # compute weights
  XtX = t(X) %*% X
  
  DXtX = diag(XtX)
  DXtX_rep = rep(DXtX, p) 
  DXtX_mat = matrix(DXtX_rep, n*p, p, byrow=FALSE)
  Diff_mat = XtX - diag(DXtX)
  
  D = matrix(1, n, n)
  for(i in 1:n){
    for(j in 1:n){
      D[i, j] = dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, .1)
    }
  }
  for(i in 1:n){
    D[, i] = n * (D[, i] / sum(D[, i])) #Scaling the weights so that they add up to n
         # D[,i]=1 # When there is no covariate information, set the weights to be 1 throughout.
  }
  
  alpha = rep(0.2, n*p) #Initialization of the inclusion probability matrix for a fixed variable, with i th row corresponding to i th subject.
  sigmabeta_sq = 3 #Initialization for hyperparameter
  mu = rep(0, p) # Variational parameter
  true_pi = 0.5 #Hyperparameter
  
  y_long_vec = as.vector(t(y %*% matrix(1, 1, p)))
  Xty=t(X) %*% as.vector(y)
  beta_mat = matrix(beta, n, p, byrow=TRUE)
  mu_mat = beta_mat
  
  D_long = matrix(0, n*p, n)
  for( i in 1:n){
    D_long[, i] = matrix(t(D[, i] %*% matrix(1, 1, p)), n*p, 1)
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
  beta_matr = matrix(0, n, p)
  
  ####################tuning hyperparameters##################################
  idmod = varbvs::varbvs(X_mat, y, Z=Z[, 1], verbose = FALSE)#Setting hyperparameter value as in Carbonetto Stephens model
  inprob = idmod$pip
  rest_index_set = setdiff(c(1:(p+1)), resp_index)
  
  sigmasq = mean(idmod$sigma)
  pi_est = mean(1 / (1 + exp(-idmod$logodds)))
  sigmavec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)
  elb1 = matrix(0, length(sigmavec), 1)
  for (j in 1:length(sigmavec)) {
    res = epwpl::cov_vsvb(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmavec[j], pi_est)
    elb1[j] = res$var.elbo
    
  }
  sigmabeta_sq = sigmavec[which.max(elb1)] #Choosing hyperparameter based on ELBO maximization
  
  result = cov_vsvb(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmabeta_sq, pi_est)
  incl_prob = result$var.alpha
  mu0_val = result$var.mu0_lambda
  
  heat_alpha = matrix(incl_prob, n, p, byrow=TRUE)
  mylist[[resp_index]] = heat_alpha
}

beta = matrix(0, p+1, p+1)
for(i in 1:(p+1)){
  for(j in 1:(p+1)){
    beta[i,j] = (Lam1[i] != 0 & Lam1[j] != 0)
  }}
diag(beta) = 0

accuracy = rbindlist(lapply(1:n, function(SUBJECT) {

  alph = matrix(0, p+1, p+1)

  for(i in 1:(p+1)) {
    alph[i, -i] = mylist[[i]][SUBJECT,]; #Individual specific inclusion probability matrix
  }
  
  a = alph
  for(i in 1:(p+1)){
    for(j in i:(p+1)){
      #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
      a[i, j] = mean(c(heat_alpha[i, j], heat_alpha[j, i]))
      a[j, i] = a[i, j]
    }
  }
  
  heat_alpha = a
  
  alphvec = sort(as.vector(heat_alpha[which(heat_alpha != 0)]))
  
  selection1 = 1*(heat_alpha > 0.5)
  
  data.table(
    sensitivity = sum(selection1 & beta) / sum(beta),
    specificity = sum(!selection1 & !beta) / sum(!beta),
    individual = SUBJECT,
    
  )
}))

data = reshape::melt(t(beta))
fig = ggplot(data, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color = "brown") +
  scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "legend" ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  labs( x=expression(bold(Variables)), y=expression(bold(Variables)), title=expression(bold("True Dependence")) ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size=15, face = "bold", colour = "black"),
        axis.title = element_text(size=30, face = "bold"),
        legend.text = element_text(face="bold",size = 25),
        legend.key.size = unit(2,'lines')) +
  coord_equal()
  

      

# Discrete covariate, dependent -------------------------------------------
n = 100
p = 10

Lam1 = c(3, 3, 3, 3, rep(0, p-3)) * 5
Lam2 = c(rep(0, p-3), 3, 3, 3, 3) * 5  #For Z[i]= 0.1 #corresponds to covariate dependent model, uncomment to try this out.

Var1 = solve(Lam1 %*% t(Lam1) + diag(rep(10, p+1))) #covariance matrix for covariate level 1
Var2 = solve(Lam2 %*% t(Lam2) + diag(rep(10, p+1))) #covariance matrix for covariate level 2

X1 = MASS::mvrnorm(n/2, rep(0, p+1), Var1)
X2 = MASS::mvrnorm(n/2, rep(0, p+1), Var2)

Z = matrix(-1, n, p) #Initializing the covariate matrix

beta = matrix(0, n, p) # Ground truth of the dependence structure
resp_index = 1;# The index we consider as response
mylist = rep(list(beta), p+1) #The variable specific inclusion probability matrix:ith row corresponds to the dependence structure for the i th subject, j th matrix corresponds to 
# the j th variable as response and the remaining as predictors.
data_mat = rbind(X1, X2)

for (i in 1:n) {
  beta[i, ] = (t(Lam1[-resp_index]) > 0)*(i<=n/2) + 
    (t(Lam2[-resp_index]) > 0)*(i>n/2) #Ground truth
  
  for(j in 1:p){
    # Z[i,j]=rnorm(1,-2,.1)*(i<50) +rnorm(1,2,0.1)*(i>=50) #Uncomment to lay with other covariate values 
    Z[i,j] = -.1*(i<=n/2)  + .1*(i>n/2)
  }
}

y = data_mat[, resp_index] #Set variable number `resp_index` as the response

X_mat = data_mat[, -resp_index] #Set the remaining p variables as predictor.
X_vec = matrix(0, n*p, 1)

X = matrix(rep(0, n^2*p), nrow=n, ncol=n*p)

for(i in 1:n){
  for(j in 1:p){
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

q = matrix(2, n, 1)

sigmasq=1 #Initialization of the hyperparameter value
E <- rnorm(n, 0, sigmasq)

XtX = t(X) %*% X

DXtX = diag(XtX)
DXtX_rep = rep(DXtX, p)
DXtX_mat = matrix(DXtX_rep, n*p, p, byrow=FALSE)
Diff_mat = XtX - diag(DXtX)

D = matrix(1, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    D[i, j] = dnorm(norm(Z[i,] - Z[j,], "2"), 0, .1)
  }
}
for (i in 1:n) {
  D[, i] = n*(D[,i] / sum(D[,i])) #Scaling the weights so that they add up to n
  #      D[,i]=1 # When there is no covariate information, set the weights to be 1 throughout.
}

