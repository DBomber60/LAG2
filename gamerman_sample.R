library(mvtnorm)

# MCMC sampling scheme for glms
# based on: https://link.springer.com/content/pdf/10.1023/A:1018509429360.pdf

# input: glm: exponential family-specific functions, current beta value
# output: mean/ variance of the N to sample from; mean = (X^T (W) X) ^1 X^T W z; variance = (X^T (W) X) ^1

getparams = function (Y, X, b, b1, b1inv, b2, beta.curr) {
  eta = X %*% beta.curr
  mu = b1(eta)
  W <- sapply(eta, b2) # b''(theta) = V(mu)
  z <- eta + (Y - mu) / W # working response
  wlm <- lm(z ~ X - 1, weights = W) # weighted least squares
  beta <- coef(wlm)
  return(list(mean = beta, sigma = summary(wlm)$cov.unscaled))
}

# input: partition matrix X; design matrix X
getparamsCOX = function (S, X, beta.curr, c) {
  p = dim(X)[2]
  e_eta = as.numeric(exp(X %*% beta.curr))
  g = sum(e_eta)
  W = sum(S)*(diag(e_eta) * g - e_eta %*% t(e_eta) )/g^2 + diag(.001, nrow = dim(X)[1])
  z = log(e_eta) + solve(W) %*% (colSums(S) - sum(S) * e_eta/g)
  sigma_mat = solve( 1/c * diag(p) + t(X) %*% W %*% X)
  mu_vector = sigma_mat %*% (t(X) %*% W %*% z)
  return(list(mean = mu_vector, sigma = sigma_mat))
}


# log-likelihood for a sample from an exponential family in which theta = eta

logpi_beta = function(Y, X, b, beta) {
  eta = X %*% beta
  return(sum(Y * eta - b(eta)))
}

# log likelihood for cox model - given partition matrix/ k-vector/ beta vector/ design matrix

logpi_surv = function(S, X, beta) {
  # S_tilde = S %*% X
  g = log(sum(X %*% beta))
  #print(g)
  return ( sum (S %*% X %*% beta - sum(S) * g ) )
}

# if survival, then Y is the partition matrix S

gamerman_mcmc = function(Y, X, beta.init, nIter = 1000, b = 1, b1 = 1, b1inv = 1, b2 = 1, surv = F, c=1) {
  
  # initialize
  beta.curr = as.numeric(beta.init)
  
  if (surv == T) {
    params.curr = getparamsCOX(Y, X, beta.curr, c)
    logpi.curr = logpi_surv(Y, X, beta.curr)
  } else {
    params.curr = getparams(Y, X, b, b1, b1inv, b2, beta.curr)
    logpi.curr = logpi_beta(Y, X, b, beta.curr)
  }
  
  # hold results in Res matrix
  p = dim(X)[2]
  Res = array(0, dim = c(nIter,p))

  for (i in 1:nIter) {
    
    # now sample a new beta
    beta.prop = as.numeric(params.curr$mean + t(chol(params.curr$sigma)) %*% rnorm(p))
    print(beta.prop)

    # get params of proposal distribution
    
    if (surv == T) {
      params.prop = getparamsCOX(Y, X, beta.prop, c)
      print(sum (X %*% beta.prop) ) # maybe just reject if sum(X %*% beta.prop) < 0 ?
      logpi.prop = logpi_surv(Y, X, beta.prop)
      #print(logpi.prop)
    } else {
      params.prop = getparams(Y, X, b, b1, b1inv, b2, beta.prop)
      logpi.prop = logpi_beta(Y, X, b, beta.prop)
    }
    
    # exp(log(q(beta^new, beta^old)) - log(q(beta^old, beta^new)))
    #qc2 = exp(dmvnorm(beta.curr, mean = params.prop$mean, sigma = params.prop$sigma, log = T) - 
    #  dmvnorm(beta.prop, mean = params.curr$mean, sigma = params.curr$sigma, log = T))
    

    acc = exp( logpi.prop - logpi.curr ) #* qc2
    
    if (runif(1) < acc) {
      #set beta to the proposed beta
      beta.curr = beta.prop
      params.curr = params.prop
      logpi.curr = logpi.prop
    }
    Res[i,] = beta.curr
  }
  return(Res)
}





