##
# Description: This file implements the spike and slab bayesian model for group feature selection 
#              
# Author: Daniel HernÃ¡ndez Lobato
#
# Year: 2011
#
#

sigmoid <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))

##
# Function: groupSpikeAndSlab
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: Trains a group slab bayesian model on the data 
#
# Arguments: 
#		X -> a matrix containing the attributes of the training data (one row per sample)
#		Y -> a vector with the class labels 1 or -1
#		tau -> Inverse of the variance of the noise
#		groups -> a vector with the group of each feature. This is a vector with postitive integers 1, 2, 3,...
#                         that identify each group.
#		p0 -> vector of prior probabilities for the gammas of each group
#		verbose -> shall we print information?
#		initializeApproximation -> approximation used to initialize EP
#
# Returns: 	A ranking of the features 
#

#' @export
GroupSpikeAndSlab <- function(X, Y, tau = 1, groups = NULL, p1 = NULL, v1 = 1, 
                              verbose = TRUE, initialApprox = NULL, opt = FALSE,
                              delta = 1e-4, max_iter = 200,
                              damping = .5, k = .99,
                              opt_lower = c(0, 0), opt_upper = c(Inf, Inf)) {
  
  # We compute some useful constants
  
  if (is.null(groups))
    groups <- 1 : ncol(X)
  
  nGroups <- max(groups)
  
  d <- ncol(X)
  
  if (is.null(p1)) 
    p0 <- rep(0.5, nGroups)
  else
    p0 <- p1
  
  # Number of samples per task
  
  n <- nrow(X)
  
  # We extend the data with one additional column to represent the biases. The prior 
  # probability of the associated gamma variable is 1 and the variance of the gaussian 100
  
  X <- t.default(X)
  
  v1 <- c(rep(v1, nGroups))
  
  # We initialize the posterior approximation
  
  posteriorApproximation <- initializeApproximation(X, Y, tau, p0, v1, groups, initialApprox)
  
  iters <- 0
  convergence <- FALSE
  
  while(iters < max_iter && ! convergence) {
    
    posteriorApproximationOld <- posteriorApproximation
    
    # We process the terms corresponding to the prior
    
    if (any(c(p0) != 1))
      posteriorApproximation <- processPriorTerms(X, Y, v1, posteriorApproximation, damping)
    
    # optimize the variance hyper parameters 
    if (opt) {
      
      mlik = function(params) {
        tau = params[1]
        v1 = params[2]
        
        value = evaluateEvidence(X, Y, p0, posteriorApproximation, tau, rep(v1, nGroups))
        return(-value)
      }
      
      hyper_opt = dfoptim::nmkb(par = c(tau, v1[1]),
                                fn = mlik,
                                lower = opt_lower,
                                upper = opt_upper)

      tau = hyper_opt$par[1]
      v1 = rep(hyper_opt$par[2], nGroups)
    }
    
    # We look for convergence in the EP algorithm
    
    convergence <- checkConvergence(posteriorApproximationOld, posteriorApproximation, 
                                    delta, verbose)
    
    iters <- iters + 1
    damping <- damping * k
  }
  
  # We evaluate the models evidence
  
  evidence <- evaluateEvidence(X, Y, p0, posteriorApproximation, tau, v1)
  
  # We compute the feature ranking 
  
  ranking <- sort(posteriorApproximation$p[ -length(posteriorApproximation$p) ], decreasing = TRUE, index.return = TRUE)$ix
  
  # We compute the mean and the variance of the marginals
  
  lambda <- posteriorApproximation$viTilde^-1
  upsilon <- posteriorApproximation$miTilde * posteriorApproximation$viTilde^-1
  
  Delta_X <- matrix(posteriorApproximation$nujTilde, d, n) * X
  vectMul <- (X %*% upsilon  + posteriorApproximation$nujTilde^-1 * posteriorApproximation$mujTilde)
  
  Inv <- solve(diag(lambda^-1) + t(X) %*% (matrix(posteriorApproximation$nujTilde, d, n) * X), tol = 1e-100)
  
  meanMarginals <- posteriorApproximation$nujTilde * vectMul - Delta_X %*% Inv %*% (t.default(Delta_X) %*% vectMul)
  varMarginals <- as.vector(posteriorApproximation$nujTilde - rowSums(Delta_X * t.default(Inv %*% t(Delta_X))))
  
  # We return the model evidence, the posterior approximation and the feature ranking
  
  list(evidence = evidence,
       posteriorApproximation = posteriorApproximation, 
       X = X, Y = Y, 
       ranking = ranking, 
       meanMarginals = meanMarginals, 
       varMarginals = varMarginals, 
       opt = opt,
       v_noise = 1 / tau, 
       v_slab = v1[1],
       iters = iters)
}

##
# Function: initializeApproximation
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: initializes the posterior approximation
#
# Arguments: 
#		X -> a matrix containing the attributes of the training data (one row per sample)
#		Y -> a matrix with the targets 
#		tau -> Inverse of the variance of the noise
#		p0 -> vector of prior probabilities for the gammas
#		v1 -> vector of variances for the slab gaussian
#		groups -> a vector with the group of each feature. This is a vector with postitive integers 
#			  1, 2, 3,... that identify each group.

# Returns: 	The posterior approximation
#

initializeApproximation <- function(X, Y, tau, p0, v1, groups, initialApprox = NULL) {
  
  
  # We set some useful constants. This is the number of tasks
  
  n <- ncol(X) 
  d <- nrow(X)
  
  # We define the paramters of the terms corresponding to the prior
  
  siTildePrior <- mujTilde <- nujTilde <- pjTilde <- rep(0, d)
  
  # We initialize the posterior approximation to the prior
  
  viTilde  <- rep(tau^-1, n)
  miTilde  <- as.vector(Y)
  siTildeLikelihood <- rep(log(1 / sqrt(2 * pi * tau^-1)), n)
  
  
  # We set the prior term for esp equal to the prior and initialize the posterior approximation
  # to the prior
  
  pjTilde <- rep(0, d)
  nujTilde <- p0[ groups ] * v1[ groups ]
  
  # Here log(1 + 1) is the normalization constant of the Bernoulli part
  
  siTildePrior <- - 0.5 * log(2 * pi * p0[ groups ] * v1[ groups ]) + log(1 + 1)
  
  # We store in a list the posterior approximation 
  
  # Some dummie variables not really used
  
  if (is.null(initialApprox))
    posteriorApproximation <- list(miTilde = miTilde, viTilde = viTilde, 
                                   siTildeLikelihood = siTildeLikelihood,
                                   siTildePrior = siTildePrior, mujTilde = mujTilde, nujTilde = nujTilde, p = logit(p0),
                                   pjTilde = pjTilde, groups = groups, p0 = logit(p0))
  else
    posteriorApproximation <- list(miTilde = miTilde, viTilde = viTilde, 
                                   siTildeLikelihood = siTildeLikelihood,
                                   siTildePrior = initialApprox$siTildePrior, 
                                   mujTilde = initialApprox$mujTilde, nujTilde = initialApprox$nujTilde, p = initialApprox$p,
                                   pjTilde = initialApprox$pjTilde, groups = groups, p0 = logit(p0))
  
  posteriorApproximation
}

##
# Function: processPriorTerms
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: process the terms of the prior
#
# Arguments: 
#		X -> a matrix containing the attributes of the training data (one row per sample)
#		Y -> a matrix with the targets 
#		v1 -> variance of the slab gaussian
#		posteriorApproximation -> posterior approximation
#		damping -> Damping factor to improve convergence. 
#
# Returns: 	The posterior approximation
#

processPriorTerms <- function(X, Y, v1, posteriorApproximation, damping) {
  
  # We set some useful constants. 
  
  n <- ncol(X)
  d <- nrow(X)
  groups <- posteriorApproximation$groups
  nGroups <- max(posteriorApproximation$groups)
  v0 <- rep(1e-10, nGroups)
  
  # We process each prior for each diffrent class label
  
  meanMarinals <- varMarginals <- rep(0,  d)
  
  # We extract the mean and the variance of the marginals
  # We compute first the lambda matrix and the upsilon vector
  
  lambda <- posteriorApproximation$viTilde^-1
  upsilon <- posteriorApproximation$miTilde * posteriorApproximation$viTilde^-1
  
  Delta_X <- matrix(posteriorApproximation$nujTilde, d, n) * X
  vectMul <- (X %*% upsilon + posteriorApproximation$nujTilde^-1 * posteriorApproximation$mujTilde)
  
  Inv <- solve(diag(lambda^-1) + t(X) %*% (matrix(posteriorApproximation$nujTilde, d, n) * X), tol = 1e-100)
  
  meanMarinals <- posteriorApproximation$nujTilde * vectMul - Delta_X %*% Inv %*% (t.default(Delta_X) %*% vectMul)
  varMarginals <- posteriorApproximation$nujTilde - colSums(((Inv) %*% t.default(Delta_X)) * t.default(Delta_X))
  
  # We remove the approximate term from the posterior approximation
  
  varMarginalsOld <- (varMarginals^-1 - posteriorApproximation$nujTilde^-1)^-1
  varMarginalsOld[ ! is.finite(varMarginalsOld) ] <- 1e6
  
  meanMarinalsOld <- varMarginalsOld * (varMarginals^-1 * meanMarinals - posteriorApproximation$nujTilde ^-1 * posteriorApproximation$mujTilde)
  
  pOld <- (posteriorApproximation$p[ groups ] - posteriorApproximation$pjTilde)  
  
  # If one term has negative variance or missing we do not perform the update
  
  indx <- which(! is.infinite(varMarginalsOld) & varMarginalsOld > 0)
  
  pOld <- pOld[ indx ]
  varMarginalsOld <- varMarginalsOld[ indx ]
  meanMarinalsOld <- meanMarinalsOld[ indx ]
  v0K <- (v0[ groups ])[ indx ]
  v1K <- (v1[ groups ])[ indx ]
  
  # Now we compute an updated posterior distribution
  
  logG1 <- dnorm(0, mean = meanMarinalsOld, sd = sqrt(varMarginalsOld + v1K), log = TRUE) 
  logG0 <- dnorm(0, mean = meanMarinalsOld, sd = sqrt(varMarginalsOld + v0K), log = TRUE)
  
  Zj <- plogis(pOld) * exp(logG1) + plogis(- pOld) * exp(logG0)
  
  c1 <- plogis(pOld + logG1 - logG0) * - meanMarinalsOld / (varMarginalsOld + v1K) + plogis(- pOld - logG1 + logG0) * 
    - meanMarinalsOld / (varMarginalsOld + v0K)
  
  c2 <- (plogis(pOld + logG1 - logG0) * (meanMarinalsOld^2 / (varMarginalsOld + v1K)^2 - 1 / (varMarginalsOld + v1K)) + 
           plogis(- pOld - logG1 + logG0) * (meanMarinalsOld^2 / (varMarginalsOld + v0K)^2 - 1 / (varMarginalsOld + v0K)))
  
  c3 <- c1^2 - c2
  
  varMarginalsNew <- varMarginalsOld - c3 * varMarginalsOld^2
  meanMarinalsNew <- meanMarinalsOld + c1 * varMarginalsOld
  
  # Now we update the approximate terms using damping factors
  
  old_nujTilde <- posteriorApproximation$nujTilde[ indx ]
  new_nujTilde <- c3^-1 - varMarginalsOld
  
  # We avoid variances equal to zero or negative
  
  new_nujTilde[ new_nujTilde == 0 ] <- 1e-10
  new_nujTilde[ new_nujTilde < 0 ] <- 1e2
  
  old_mujTilde <- posteriorApproximation$mujTilde[ indx ]
  new_mujTilde <- meanMarinalsOld + c1 * (new_nujTilde + varMarginalsOld)
  old_pjTilde <- posteriorApproximation$pjTilde[ indx ]
  new_pjTilde <- logG1 - logG0
  
  posteriorApproximation$nujTilde[ indx ] <- (damping * new_nujTilde^-1 + (1 - damping) * old_nujTilde^-1)^-1  
  posteriorApproximation$mujTilde[ indx ] <- posteriorApproximation$nujTilde[ indx ] * 
    (damping * new_nujTilde^-1 * new_mujTilde + (1 - damping) * old_nujTilde^-1 * old_mujTilde)
  posteriorApproximation$pjTilde[ indx ] <- damping * new_pjTilde + (1 - damping) * old_pjTilde
  
  posteriorApproximation$siTildePrior[ indx ] <- log(exp(logG0) + exp(logG1)) + 0.5 * 
    log(1 + varMarginalsOld * posteriorApproximation$nujTilde[ indx ]^-1) + (0.5 * c1^2 / c3)
  
  # We update the posterior distribution. We take care of non-updated groups
  
  iGroup <- sort(unique(groups[ indx ]))
  
  # We avoid overflows while recomputing the posterior approximation as the normalized product of the ~ti
  
  posteriorApproximation$p[ iGroup ] <- tapply(posteriorApproximation$pjTilde[ indx ], 
                                               as.factor(posteriorApproximation$groups[ indx ]), sum) + posteriorApproximation$p0[ iGroup ]
  
  # Now we update the matrix A and B and the vecotr h
  
  posteriorApproximation
}

##
# Function: checkConvergence
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: looks for convergence in the EP algorithm
#
# Arguments: 
#		posteriorApproximationOld -> posterior approximation
#		posteriorApproximationNew -> posterior approximation
#		verbose -> shall we print information?
#
# Returns: 	A boolean indicating convergence
#

checkConvergence <- function(posteriorApproximationOld, posteriorApproximationNew, 
                             delta, verbose) {
  
  # We evaluate the maximum change within the posterior approximation
  
  maxChange <- 0
  
  maxChange <- max(maxChange, abs(posteriorApproximationOld$A - posteriorApproximationNew$A))
  maxChange <- max(maxChange, abs(posteriorApproximationOld$h - posteriorApproximationNew$h))
  maxChange <- max(maxChange, abs(posteriorApproximationOld$p - posteriorApproximationNew$p))
  
  if (verbose)
    cat("EP: max change", maxChange, "\n")
  
  if (maxChange < delta)
    TRUE
  else
    FALSE
}

##
# Function: evaluateEvidence
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: evaluates the models evidence provided by the EP algorithm
#
# Arguments: 
#		X -> a matrix containing the attributes of the training data (one row per sample)
#		Y -> a matrix with the targets
#		p0 -> vector of prior probabilities for the gammas
#		posteriorApproximation -> posterior approximation
#
# Returns: 	Log evidence of the model
#

evaluateEvidence <- function(X, Y, p0, posteriorApproximation, tau, v1) {
  
  # We set some useful constants. This is the number of tasks
  
  
  n <- ncol(X) 
  d <- nrow(X)
  
  # We compute the log evidence step by step
  
  prodActive <- tapply(plogis(posteriorApproximation$pjTilde), 
                       as.factor(posteriorApproximation$groups), prod) * plogis(posteriorApproximation$p0)
  prodInactive <- tapply(plogis(- posteriorApproximation$pjTilde), 
                         as.factor(posteriorApproximation$groups), prod) * plogis(- posteriorApproximation$p0)
  
  D_j <- prodActive + prodInactive
  
  C <- (diag(posteriorApproximation$viTilde) + t.default(X) %*% (matrix(posteriorApproximation$nujTilde, d, n) * X))
  
  gaussianPartOriginal <- -0.5 * (n * log(2 * pi) + determinant(C)$modulus[[ 1 ]] + 
                                    (t.default(Y - t(X) %*% posteriorApproximation$mujTilde) %*% solve(C) %*% (Y - t(X) %*% posteriorApproximation$mujTilde))[ 1, 1 ]) - 
    sum(-0.5 * log(2 * pi * posteriorApproximation$nujTilde))
  
  # We add the contribution from the beta part of each approximate ti
  
  log_s1 <- - n / 2 * log(2 * pi * tau^-1) - 1 / 2 * tau * sum(Y * Y)
  
  n <- ncol(X)
  d <- nrow(X)
  groups <- posteriorApproximation$groups
  nGroups <- max(posteriorApproximation$groups)
  v0 <- rep(1e-10, nGroups)
  
  # We process each prior for each diffrent class label
  
  meanMarinals <- varMarginals <- rep(0,  d)
  
  # We extract the mean and the variance of the marginals
  # We compute first the lambda matrix and the upsilon vector
  
  lambda <- posteriorApproximation$viTilde^-1
  upsilon <- posteriorApproximation$miTilde * posteriorApproximation$viTilde^-1
  
  Delta_X <- matrix(posteriorApproximation$nujTilde, d, n) * X
  vectMul <- (X %*% upsilon + posteriorApproximation$nujTilde^-1 * posteriorApproximation$mujTilde)
  
  Inv <- solve(diag(lambda^-1) + t(X) %*% (matrix(posteriorApproximation$nujTilde, d, n) * X), tol = 1e-100)
  
  meanMarinals <- posteriorApproximation$nujTilde * vectMul - Delta_X %*% Inv %*% (t.default(Delta_X) %*% vectMul)
  varMarginals <- posteriorApproximation$nujTilde - colSums(((Inv) %*% t.default(Delta_X)) * t.default(Delta_X))
  
  # We remove the approximate term from the posterior approximation
  
  varMarginalsOld <- (varMarginals^-1 - posteriorApproximation$nujTilde^-1)^-1
  varMarginalsOld[ ! is.finite(varMarginalsOld) ] <- 1e6
  
  meanMarinalsOld <- varMarginalsOld * (varMarginals^-1 * meanMarinals - posteriorApproximation$nujTilde ^-1 * posteriorApproximation$mujTilde)
  
  pOld <- (posteriorApproximation$p[ groups ] - posteriorApproximation$pjTilde)  
  
  # If one term has negative variance or missing we do not perform the update
  
  indx <- which(! is.infinite(varMarginalsOld) & varMarginalsOld > 0)
  
  pOld <- pOld[ indx ]
  varMarginalsOld <- varMarginalsOld[ indx ]
  meanMarinalsOld <- meanMarinalsOld[ indx ]
  v0K <- (v0[ groups ])[ indx ]
  v1K <- (v1[ groups ])[ indx ]
  
  # Now we compute an updated posterior distribution
  
  logG1 <- dnorm(0, mean = meanMarinalsOld, sd = sqrt(varMarginalsOld + v1K), log = TRUE) 
  logG0 <- dnorm(0, mean = meanMarinalsOld, sd = sqrt(varMarginalsOld + v0K), log = TRUE)
  
  Zj <- plogis(pOld) * exp(logG1) + plogis(- pOld) * exp(logG0)
  
  c1 <- plogis(pOld + logG1 - logG0) * - meanMarinalsOld / (varMarginalsOld + v1K) + 
    plogis(- pOld - logG1 + logG0) * - meanMarinalsOld / (varMarginalsOld + v0K)
  
  c2 <- (plogis(pOld + logG1 - logG0) * (meanMarinalsOld^2 / (varMarginalsOld + v1K)^2 - 1 / (varMarginalsOld + v1K)) + 
           plogis(- pOld - logG1 + logG0) * (meanMarinalsOld^2 / (varMarginalsOld + v0K)^2 - 1 / (varMarginalsOld + v0K)))
  
  c3 <- c1^2 - c2
  
  log_s2j <- log(Zj) - log(plogis(pOld) * plogis(posteriorApproximation$pjTilde) + 
                             plogis(-pOld) * plogis(-posteriorApproximation$pjTilde)) - 0.5 * log(2 * pi * posteriorApproximation$nujTilde) - 
    dnorm(0, meanMarinalsOld - posteriorApproximation$mujTilde, 
          sqrt(posteriorApproximation$nujTilde + varMarginalsOld), log = TRUE)
  
  log_s2 <- sum(log_s2j)
  
  log_s3 <- 0
  
  prodActive <- tapply(plogis(posteriorApproximation$pjTilde), 
                       as.factor(posteriorApproximation$groups), prod) * plogis(posteriorApproximation$p0)
  prodInactive <- tapply(plogis(- posteriorApproximation$pjTilde), 
                         as.factor(posteriorApproximation$groups), prod) * plogis(- posteriorApproximation$p0)
  
  D_j <- prodActive + prodInactive
  
  upsilon <- tau * X %*% Y + posteriorApproximation$mujTilde * posteriorApproximation$nujTilde^-1
  
  Xv <- matrix(upsilon, nrow = 1)
  
  Xs <- t(X) * matrix(posteriorApproximation$nujTilde, n, d, byrow = TRUE)
  M <- solve(diag(n) * tau^-1 + Xs %*% X)
  
  Tv <- Xs %*% t.default(Xv)
  Cv <- t.default(Xv) * matrix(posteriorApproximation$nujTilde, d, nrow(Xv)) - t(Xs) %*% (M %*% Tv)
  values <- colSums(Cv * t.default(Xv))
  
  gaussianPart <- log_s1 + 
    - 1 / 2 * sum(posteriorApproximation$mujTilde^2 * posteriorApproximation$nujTilde^-1) +
    + d / 2 * log(2 * pi) + 1 / 2 * values + 
    + 1 / 2 * (sum(log(posteriorApproximation$nujTilde)) - determinant(diag(n) + tau *
                                                                         t.default(X) %*% (matrix(posteriorApproximation$nujTilde, d, n) * X))$modulus[[ 1 ]])
  
  evidence <- log_s2 + sum(log(D_j)) + gaussianPart 
  evidence2 <- log_s2 + sum(log(D_j)) + gaussianPartOriginal
  
  evidence
}

##
# Function: predictSpikeAndSlab
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: computes the prediction for new instances
#
# Arguments: 
#		EPApproximation -> posterior approximation
#		x -> a matrix containing the attributes of the test data (one column per sample)
#
# Returns: 	Class labels of the different instances
#

predictSpikeAndSlab <- function(EPApproximation, x) {
  
  # We check if we have more than one instance
  
  if (!is.matrix(x))
    x <- matrix(x, 1, length(x))
  
  x <- t.default(x)
  
  # We compute some useful constants
  
  n <- ncol(EPApproximation$X) 
  d <- nrow(EPApproximation$X)
  
  nSamplesToPredict  <- ncol(x)
  
  lambda <- EPApproximation$posteriorApproximation$viTilde^-1
  upsilon <- EPApproximation$posteriorApproximation$miTilde * EPApproximation$posteriorApproximation$viTilde^-1
  
  Delta_X <- matrix(EPApproximation$posteriorApproximation$nujTilde, d, n) * EPApproximation$X
  vectMul <- (EPApproximation$X %*% upsilon  + EPApproximation$posteriorApproximation$nujTilde^-1 * EPApproximation$posteriorApproximation$mujTilde)
  
  Inv <- solve(diag(lambda^-1) + t(EPApproximation$X) %*% (matrix(EPApproximation$posteriorApproximation$nujTilde, d, n) * EPApproximation$X), tol = 1e-100)
  
  meanMarginals <- EPApproximation$posteriorApproximation$nujTilde * vectMul - Delta_X %*% Inv %*% (t.default(Delta_X) %*% vectMul)
  varMarginals <- as.vector(EPApproximation$posteriorApproximation$nujTilde - rowSums(Delta_X * t.default(Inv %*% t(Delta_X))))
  
  mkx <- t.default(x) %*% meanMarginals
  
  Delta_Xnew <- t.default(x) %*% (matrix(EPApproximation$posteriorApproximation$nujTilde, d, n) * EPApproximation$X)
  xVkx <- colSums(matrix(EPApproximation$posteriorApproximation$nujTilde, d, nSamplesToPredict) * x^2) - colSums(((Inv) %*% t.default(Delta_Xnew)) * t.default(Delta_Xnew))
  
  list(mean = mkx, var = xVkx)
}

##
# Function: estimatePredictivePerformance
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description: computes the prediction for new instances
#
# Arguments: 
#		EPApproximation -> posterior approximation, as returned by the multiClassSpikeAndSlabFeatureRanking
#		x -> a matrix containing the attributes of the training data (one column per sample)
#
# Returns: 	Class labels of the different instances
#

estimatePredictivePerformance <- function(EPApproximation, sdY = 1, meanY = 0, nFolds = 10) {
  
  x <- EPApproximation$X
  
  e <- rep(0, 10)
  
  # We compute some useful constants
  
  n <- ncol(EPApproximation$X) 
  d <- nrow(EPApproximation$X)
  
  nSamplesToPredict  <- ncol(x)
  
  CV <- list()
  
  perm <- 1 : n
  tablas <- tapply(1 : n, cut(1 : n, breaks = nFolds), function(x) perm[ x ])
  CV <- tablas
  
  for (i in 1 : nFolds) {
    
    # We set to uniform the first term
    
    lambda <- EPApproximation$posteriorApproximation$viTilde[ -CV[[ i ]] ]^-1
    upsilon <- EPApproximation$posteriorApproximation$miTilde[ -CV[[ i ]] ] * EPApproximation$posteriorApproximation$viTilde[ -CV[[ i ]] ]^-1
    
    Delta_X <- matrix(EPApproximation$posteriorApproximation$nujTilde, d, n - length(CV[[ i ]])) * EPApproximation$X[, -CV[[ i ]] ]
    vectMul <- (EPApproximation$X[ , -CV[[ i ]] ] %*% upsilon + EPApproximation$posteriorApproximation$nujTilde^-1 * EPApproximation$posteriorApproximation$mujTilde)
    
    Inv <- solve(diag(lambda^-1) + t(EPApproximation$X[, -CV[[ i ]] ]) %*% (matrix(EPApproximation$posteriorApproximation$nujTilde, d, n - length(CV[[ i ]])) * 
                                                                              EPApproximation$X[ , -CV[[ i ]] ]), tol = 1e-100)
    
    meanMarginals <- EPApproximation$posteriorApproximation$nujTilde * vectMul - Delta_X %*% Inv %*% (t.default(Delta_X) %*% vectMul)
    
    e[ i ] <- mean(((EPApproximation$Y[ CV[[ i ]] ] * sdY + meanY) - (t(x[, CV[[ i ]] ]) %*% meanMarginals* sdY + meanY))^2)
  }
  
  mean(e)
}




##
# Function: tuneParam
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description:  finds the best hyperparameters for the  model
#
# Arguments: 
#               X <- matrix with the different data 
#               Y <- matrix with the different labels 
#               groupIndex <- groups of each attribute
#               nFolds <- number of folds used in the validation
#               verbose <- Should additional output be printed?
#
# Returns: optimal hyperparameters and error estimate
#

tuneParam <- function(X, Y, tau = 1, groupIndex, nFolds = 10, p1 = c(0.1, 0.15, 0.2, 0.25, 0.3), 
                      v1 = c(0.025, 0.05, 0.075, 0.1, 0.125), verbose = FALSE, fast = FALSE) {
  
  n <- nrow(X)	
  
  # We partition each data set using 10-fold cross validatoin
  
  CV <- list()
  
  perm <- sample(1 : n)
  tablas <- tapply(1 : n, cut(1 : n, breaks = nFolds), function(x) perm[-x])
  CV <- tablas
  
  if (verbose == TRUE)
    cat("Generating sequence.\n")
  
  
  errors <- matrix(0, length(p1), length(v1))
  
  # We do the performance assessment
  
  
  if (! fast) {
    
    for (i in 1 : nFolds) {
      
      XTrainFold <- X[ CV[[ i ]], ]
      YTrainFold <- Y[ CV[[ i ]] ]
      
      XTestFold <- X[ -CV[[ i ]], ]
      YTestFold <- Y[ -CV[[ i ]] ]
      
      if (verbose == TRUE)
        cat("Fold:", i, "Tunning for params\n")
      
      ret <- normTrain(XTrainFold, robust = FALSE)
      XTrainFold <- ret$values
      XTestFold <- normTest(XTestFold, ret$ms, ret$sds)
      
      ret <- normTrain(YTrainFold, robust = FALSE)
      YTrainFold <- ret$values
      
      for (p in 1 : length(p1))  {
        
        for (v in 1 : length(v1)) {
          
          cat("Testing p1:", p1[ p ], "v1:", v1[ v ], "\n")
          
          model <- GroupSpikeAndSlab(XTrainFold, YTrainFold, tau, groupIndex, rep(p1[ p ], 
                                                                                  length(unique(groupIndex))), v1[ v ], verbose = FALSE)
          
          # We do the test validation. 
          
          errors[ p, v ] <- errors[ p, v ] +  mean((YTestFold - (XTestFold %*% model$meanMarginals * ret$sds + ret$ms))^2)
        }
      }
    }
    
    errors <- errors / nFolds
    
  } else {
    
    XTrainFold <- X
    YTrainFold <- Y
    
    ret <- normTrain(XTrainFold, robust = FALSE)
    XTrainFold <- ret$values
    
    ret <- normTrain(YTrainFold, robust = FALSE)
    YTrainFold <- ret$values
    
    for (p in 1 : length(p1))  {
      
      for (v in 1 : length(v1)) {
        
        cat("Testing p1:", p1[ p ], "v1:", v1[ v ], "\n")
        
        model <- GroupSpikeAndSlab(XTrainFold, YTrainFold, tau, groupIndex, rep(p1[ p ], length(unique(groupIndex))), v1[ v ], verbose = FALSE)
        
        # We do the test validation
        
        errors[ p, v ] <- estimatePredictivePerformance(model, ret$sds, ret$ms)
      }
    }
  }
  
  
  res <- which.min(c(errors))
  
  v1Best <- as.integer((res - 1) / length(p1)) + 1
  p1Best <- (res - 1) %% length(p1) + 1
  
  # We compute the error for each hyperparameter 
  
  list(v1 = v1[ v1Best ], p1 = p1[ p1Best ], value = c(errors)[ res ])
}


##
# Function: tuneParamEvidence
#
# Author: Daniel HernÃ¡ndez-Lobato
#
# Description:  finds the best hyperparameter for the  model using the evidence
#
# Arguments: 
#               X <- matrix with the different data 
#               Y <- matrix with the targets 
#               groupIndex <- groups of each attribute
#               nFolds <- number of folds used in the validation
#               verbose <- Should additional output be printed?
#
# Returns: optimal hyperparameter and error estimate
#

tuneParamEvidence <- function(X, Y, tau = 1, groupIndex, p1 = c(0.1, 0.15, 0.2, 0.25, 0.3), 
                              v1 = c(0.025, 0.05, 0.075, 0.1, 0.125), verbose = FALSE, fast = FALSE) {
  
  n <- nrow(X)	
  
  XTrainFold <- X
  YTrainFold <- Y
  
  # We partition each data set using 10-fold cross validatoin
  
  if (verbose == TRUE)
    cat("Generating sequence.\n")
  
  evidences <- matrix(0, length(p1), length(v1))
  
  # We do the performance assessment
  
  for (p in 1 : length(p1))  {
    
    for (v in 1 : length(v1)) {
      
      cat("Testing p1:", p1[ p ], "v1:", v1[ v ], "\n")
      
      model <- GroupSpikeAndSlab(XTrainFold, YTrainFold, tau, groupIndex, rep(p1[ p ], 
                                                                              length(unique(groupIndex))), v1[ v ], verbose = FALSE)
      
      # We do the test validation.
      
      evidences[ p, v ] <- model$evidence
    }
  }
  
  res <- which.max(c(evidences))
  
  v1Best <- as.integer((res - 1) / length(p1)) + 1
  p1Best <- (res - 1) %% length(p1) + 1
  
  list(v1 = v1[ v1Best ], p1 = p1[ p1Best ], value = c(evidences)[ res ])
}


# These functions are used to normalize the data. The roubst parameter is not used in practice.

normTrain <- function(x, robust = FALSE){
  
  if (is.vector(x))
    x <- matrix(x, nrow = length(x), ncol = 1)
  
  means <- apply(x,2,mean)
  sds <- apply(x, 2, sd)
  
  if (any(sds == 0)) {
    means[ sds == 0 ] <- 0
    sds[ sds == 0 ] <- 1
  }
  
  ms <- means 
  sds <- sds
  
  x <- scale(x, ms, sds)
  
  list(ms=ms,sds=sds,values=x)
}

normTest<-function(x,ms,sds){ scale(x, ms, sds) }
