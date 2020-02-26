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

# log-likelihood for a sample from an exponential family in which theta = eta

logpi_beta = function(Y, X, b, beta) {
  eta = X %*% beta
  return(sum(Y * eta - b(eta)))
}


gamerman_mcmc = function(Y, X, beta.init, nIter, b, b1, b1inv, b2) {
  
  # initialize
  beta.curr = beta.init
  params.curr = getparams(Y, X, b, b1, b1inv, b2, beta.curr)
  logpi.curr = logpi_beta(Y, X, b, beta.curr)
  
  # hold results in Res matrix
  p = dim(X)[2]
  Res = array(0, dim = c(nIter,p))

  for (i in 1:nIter) {
    
    # now sample a new beta
    beta.prop = as.numeric(params.curr$mean + t(chol(params.curr$sigma)) %*% rnorm(p))
    params.prop = getparams(Y, X, b, b1, b1inv, b2, beta.prop)
    
    # log(p(beta^new))
    logpi.prop = logpi_beta(Y, X, b, beta.prop)
    
    # exp(log(q(beta^new, beta^old)) - log(q(beta^old, beta^new)))
    qc2 = exp(dmvnorm(beta.curr, mean = params.prop$mean, sigma = params.prop$sigma, log = T) - 
      dmvnorm(beta.prop, mean = params.curr$mean, sigma = params.curr$sigma, log = T))

    acc = exp( logpi.prop - logpi.curr ) * qc
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





