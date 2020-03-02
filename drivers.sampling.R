library(igraph)

############### LAG_SAMPLE: Obtain sample from model ####################

# input - design matrix (X), number of cliques to sample (k), cardinality coefficients (alpha), popularity coefficients (beta)
# output - single transaction - non-overlapping cliques sample with weight proportional to X %*% (alpha,beta)
sample_transaction = function(X, k, alpha, beta, cn) { 
  s = array(0, dim = c(1,nrow(X))) # partition vector - one column for each clique in G
  rsk = exp(X %*% c(alpha, beta))
  sampled = sample.int(length(rsk), size = k, prob = rsk)
  s[sampled] = 1
  return(list(s = s, transaction = s %*% X))
}

# sample an association graph 
# input: n, gamma
# output: G, clique design matrix based on sampled G
g_sample = function(n, gamma) {
  edges <- c()
  for (i in 1:(n - 1)) { #    for each
    for (j in (i + 1):n) { # vertex pair
      if (rbinom(1, 1, ilogit(gamma[i] + gamma[j])) == 1)
        edges <- c(edges, i, j)
    }
  }
  g1 <- make_empty_graph(n, directed = FALSE) %>% add_edges(edges)
  chordalg1 = is_chordal(g1, fillin = T, newgraph = T)$newgraph
  #g3 = chordalg2$newgraph # ghat
  return(chordalg1)
}

# make design matrix (X) given a clique set, clique number, alpha, beta
makeDesign = function(C, cn, alpha, beta) {
  # build design matrix (based on the clique set, C, of the association graph)
  X <- matrix(0, nrow = length(C), ncol = cn + n) # design
  for (i in seq_along(C)) {
    kc <- length(C[[i]]) - 1
    if (kc > 0) X[i, kc] <- 1
    X[i, cn + C[[i]]] <- 1 / (kc + 1)
  }
  return(X)
}

# input: G, k, X, cn, nt, theta, beta (n tuple), alpha (depends on graph sample)
# output: S, D

lag_sample = function(G, k, X, cn, nt, theta, gamma, beta, alpha) {
  n = gorder(G)
  S_tilde = array(0, dim = c(nt,n+cn)) # S %*% X
  S = array(0, dim = c(nt, nrow(X))) # partition matrix
  # sample 'nt' transactions
  for (i in 1:nt) {
    st = sample_transaction(X, k[i], alpha, beta, cn)
    S_tilde[i,] = st$transaction
    S[i,] = st$s
  }
  lhood.true = sum(S_tilde %*% c(alpha,beta) - k * log(sum(exp(X %*% c(alpha,beta)))) + k * log(theta))
  D = (t(apply(S_tilde,1,function(x) ifelse(x>0,1,0))))[,cn+1:n] # observed data
  return(list(D=D, lhood.true=lhood.true, k = k, S=S, S_tilde = S_tilde))
}


