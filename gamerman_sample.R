#setwd("~/Documents/LAG Model/bayesglm")
#source('irls.r')

set.seed(1)
# generate some data
n = 200
p = 2
X = matrix(rnorm(n*2), nrow=n, ncol=p)
beta.true = c(1,2)
Y = rpois(n, lambda = exp(X %*% beta.true))


#poisson.irls(Y~X-1)

# compare mcmc with frequentist
# one iteration

getparams = function (Y, X, b, b1, b1inv, b2, beta.curr) {
  eta = X %*% beta.curr
  mu = b1(eta)
  W <- sapply(eta, b2) # b''(theta) = V(mu)
  z <- eta + (Y - mu) / W # working response
  wlm <- lm(z ~ X - 1, weights = W) # weighted least squares
  beta <- coef(wlm)
  return(list(mean = beta, sigma = summary(wlm)$cov.unscaled))
}

# log-likelihood for a given beta value - poisson likelihood
logpi_beta = function(Y, X, beta) {
  eta = X %*% beta
  return(sum(Y * eta - exp(eta)))
}

qcomp = function(new_params, old_params, beta.prop, beta.curr) {
  num = -.5 * log(det(new_params$sigma)) - .5 * t(beta.curr-new_params$mean) %*% solve(new_params$sigma) %*% (beta.curr-new_params$mean)
  den = -.5 * log(det(old_params$sigma)) - .5 * t(beta.prop-old_params$mean) %*% solve(old_params$sigma) %*% (beta.prop-old_params$mean)
  return(exp(num-den))
}

# initialize
beta.curr = glm(Y~X-1, family = poisson)$coef
params.curr = getparams(Y, X, exp, exp, log, exp, beta.curr)
logpi.curr = logpi_beta(Y, X, beta.curr)

nIter = 1000

Res = array(0, dim = c(nIter,2))

for (i in 1:nIter) {
  
  # now sample a new beta
  beta.prop = params.curr$mean + t(chol(params.curr$sigma)) %*% rnorm(p)
  params.prop = getparams(Y, X, exp, exp, log, exp, beta.prop)
  
  # accept it?
  logpi.prop = logpi_beta(Y, X, beta.prop)
  #print(logpi.prop)
  
  # exp( log(q(x*,x)) - log(q(x,x*)) )
  qc = qcomp(params.prop, params.curr, beta.prop, beta.curr)
  
  acc = exp( logpi.prop - logpi.curr ) * qc
  if (runif(1) < acc) {
    #set beta to the proposed beta
    print("accept")
    beta.curr = beta.prop
    params.curr = params.prop
  }
  beta.curr = beta.prop
  params.curr = params.prop
  Res[i,] = beta.curr
}



apply(Res,2,sd)
apply(Res,2,mean)
summary(glm(Y ~ X -1, family = poisson))




