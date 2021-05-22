
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

for (i in 1:n) {
  beta[i, ] = (t(Lam1[-resp_index]) > 0)*(i<=n/2) + 
    (t(Lam2[-resp_index]) > 0)*(i>n/2) #Ground truth
  
  for(j in 1:p){
    # Z[i,j]=rnorm(1,-2,.1)*(i<50) +rnorm(1,2,0.1)*(i>=50) #Uncomment to lay with other covariate values 
    Z[i,j] = -.1*(i<=n/2)  + .1*(i>n/2)
  }
}

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

