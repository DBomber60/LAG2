source('gamerman_sample.R')
source('drivers.sampling.R')
source('rewens.R')


# test gamerman MCMC algo
# 1. poisson - good 

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
summary(glm(Y~X-1, family = poisson))


# 2. Ewens - good 

source('rewens.R')

n = 200
Y = rewens(n, theta = rep(0.25,n))
X = array(1, dim = c(n,1))

# EWENS specific canonical parameter/ cumulant function/ variance
b  = function(theta) {sum(log(exp(theta)+seq(0,n-1)))}
b1 = function(theta) {sum(exp(theta)/(seq(0,n-1)+exp(theta)))}
b2 = function(theta) {sum(seq(1,n-1)*exp(theta)/((seq(1,n-1)+exp(theta))^2))}


R = gamerman_mcmc(Y, X, beta.init = 1, nIter = 1000, b = b, b1 = b1, b1inv = log, b2 = b2)
apply(R,2,sd)
exp(apply(R,2,mean))

# ok, now for the survival part ...

# first, let's sample from the model

set.seed(1)
n = 30 # items
nt = 500 # transactions
gamma <- rnorm(n, -2) # graph vertex coefficients
theta=.2

g <- g_sample(n, gamma)
C <- lapply(cliques(g), as.vector) # quick check: `table(sapply(C, length))`
cn <- max(sapply(C, length)) - 1 # clique number of `g` (minus 1)

k = rewens(n, rep(theta,nt))

# parameters
alpha <- (1:cn) # cardinality coefficients 
beta_true <- rnorm(n) # item (vertex) coefficients
X = makeDesign(C, cn, alpha, beta = beta_true)

sampled = lag_sample(G=g, k, X, cn=cn, nt=nt, theta=theta, gamma = gamma, beta = beta_true, alpha = alpha) 

# estimate initial graph
M = t(sampled$D) %*% sampled$D
ghat = g_estimate(sampled$D, M=M, nt = nt, pval = .001)
kshat = kestimate(ghat$graph, sampled$D, cn = ghat$cn)

# let's estimate alpha, beta assuming we know S/ k













