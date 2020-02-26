source('gamerman_sample.R')


# test gamerman MCMC algo

set.seed(1)
n = 200
p = 2
X = matrix(rnorm(n*p), nrow=n, ncol=p)
beta.true = c(1,2)
Y = rpois(n, lambda = exp(X %*% beta.true))
beta.curr = glm(Y~X-1, family = poisson)$coef
R = gamerman_mcmc(Y, X, beta.curr, nIter = 1000, b = exp, b1 = exp, b1inv = log, b2 = exp)

apply(Res,2,sd)
apply(Res,2,mean)

# test with the expo families used in our model