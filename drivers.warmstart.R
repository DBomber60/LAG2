#library(arm)
#library(rlist)

#ilogit <- function (x) 1 / (1 + exp(-x)) # plogis

########## EWENS GLM ##########

irls <- function (f, b, b1, b1inv, b2, init, tol = 1e-7) {
  # initialize
  y <- model.response(model.frame(f))
  X <- model.matrix(f) # design matrix
  mu <- init(y)
  eta <- b1inv(mu) # g(mu)
  lhood <- sum(y * eta - sapply(eta,b))
  
  # iterate
  repeat {
    W <- sapply(eta, b2) # b''(theta) = V(mu)
    z <- eta + (y - mu) / W # working response
    wlm <- lm(z ~ X - 1, weights = W) # weighted least squares
    beta <- coef(wlm)
    eta <- fitted(wlm) # X %*% beta
    lhood.new <- sum(y * eta - sapply(eta,b))
    if (abs((lhood.new - lhood) / lhood) < tol) break # converged?
    lhood <- lhood.new
    mu <- sapply(eta, b1) # b'(theta)
  }
  
  # report
  list(coef = beta, var = summary(wlm)$cov.unscaled, lhood=lhood)
}

# EWENS specific canonical parameter/ cumulant function/ variance
b  = function(theta) {sum(log(exp(theta)+seq(0,n-1)))}
b1 = function(theta) {sum(exp(theta)/(seq(0,n-1)+exp(theta)))}
b2 = function(theta) {sum(seq(1,n-1)*exp(theta)/((seq(1,n-1)+exp(theta))^2))}

ewens.irls <- function (f, tol = 1e-6)
  irls(f, b=b, b1=b1, b1inv=log, b2 = b2, function (y) y + 0.5, tol)


######################## NR FITTING PROCEDURE ################################
gf = function(beta, X) {sum(exp(X %*% beta))}
lhood = function(S,beta,k,X) return(sum(S %*% beta - k * log(gf(beta, X))))

# reduced form
NRfit2 = function(S, X1, k, tol = 1e-6) {
  beta = rep(0,dim(X1)[2]) # initialize beta
  lhoodarray = lhood(S,beta,k,X1)
  loglik = lhood(S,beta,k,X1)

  repeat { 
    eta = as.numeric(exp(X1 %*% beta))
    g = sum(eta)
    W = (t(X) %*% (diag(eta) * g - eta %*% t(eta) ) %*% X)/g^2 + diag(.001, nrow = dim(X1)[2])
    U = chol(W)
    Uinv = solve(U) 
    Winv = Uinv %*% t(Uinv)
    beta_tilde = beta + sum(k)^-1 * Winv %*% (colSums(S)-sum(k)* (t(X) %*% eta)/g)
    beta_tilde[1:cn] = 1:cn # (dim(X1)[2]-n)
    
    # does this estimate increase the likelihood?
    if (lhood(S, beta_tilde, k, X1) > loglik) {beta = beta_tilde
    } else {beta = (beta*.8+beta_tilde*.2)}
    
    lhood.new = lhood(S,beta,k,X1) #sum(S %*% beta_tilde - k * log(g(beta_tilde, X1))) # lhood
    
    if (abs((lhood.new-loglik)/loglik)<tol) break # converged?
    loglik = lhood.new  
    lhoodarray = c(lhoodarray,loglik)
  }
  list(coef = beta, se = (sqrt(diag(1/sum(k) * Winv))), lhood = loglik)
}



################## COMPUTE CLIQUE COVER NUMBER/ COMPONENT CLIQUES ##################

# clique cover function
# input: chordal graph WITH NAMED VERTICES
# output: list whose elements are cliques that compose a minimum clique cover of the input graph
clique_cover = function(g) {
  v = V(g)
  peo = v[max_cardinality(g)$alpham1] # perfect elimination ordering (numeric)
  X = list()
  y = list() # stable set
  Y=list() # list of cliques (Y_i = {y_i} U {X_i})
  v1 = peo[1]
  y[[1]] = as.numeric(labels(v1))   #(labels(v1)) #as.numeric(labels(v1))
  X[[1]] = as.numeric(labels(neighbors(g,v1)))       #(labels(neighbors(g,v1)))
  Y[[1]] = c(y[[1]], X[[1]])
  
  # if all the vertices are in the first maximal clique, clique cover number is 1, don't proceed
  if(vcount(g)>length(Y[[1]])) {
    for(i in 2:(length(peo))) {
      if (!(as.numeric(labels(peo[i])) %in% unlist(X))) { #labels
        v_i = as.numeric(labels(peo[i])) # labels
        x_i = setdiff(as.numeric(labels(neighbors(g,peo[i]))), unlist(Y)) #as.numeric(labels((peo[1:(i-1)]))) ) #labels
        #x_i = as.numeric(labels(neighbors(g,peo[i]))) #as.numeric(labels((peo[1:(i-1)]))) )
        y = list.append(y, v_i)
        X = list.append(X, x_i)
        Y = list.append(Y, c(v_i,x_i))
      } 
    }
  }
  return(Y)
}

# input: clique cover of a transaction
# corresponding row in S matrix

# currently - if the clique size exceeds the current value of cn, we don't put in an alpha coefficient - FIX THIS

make_S_row = function(cc, cn) {
  S_row = array(0, dim = c(1,n+cn))
  for (j in 1:length(cc)) {
    alpha = length(cc[[j]]) - 1
    if (alpha > 0 & alpha <= cn ) { S_row[alpha] =  S_row[alpha] + 1 }
    S_row[ cc[[j]] + cn ] = S_row[ cc[[j]] + cn ] + 1/(alpha + 1)
  }
  return(S_row)
}

# input: chordal graph, data
# output: S matrix and 'k' estimate
kestimate = function(g, D, cn) {
  nt = ifelse(is.null(dim(D)[1]),1,dim(D)[1])
  S = array(0, dim = c(nt,n+cn))
  #C_list = list()
  k_est = c()
  
  if (nt == 1) {D = matrix(D,ncol = n)}
  
  for(i in 1:nt) {
    v_set = which(D[i,]>0) # vertices in first transaction (what if D is a single row ????)
    g_prime = induced.subgraph(g, v_set)
    V(g_prime)$name = v_set
    cc = clique_cover(g_prime)
    k_est[i] = length(cc)
    S[i,] = make_S_row(cc, cn)
  }
  return(list(S=S, k=k_est))
}

# input: edge, n=number of vertices
# output: associated row number of design matrix/ response vector for gamma hat estimation
edge_index = function(edge, n) {
  if(edge[1] == 1) {return(edge[2]-1)
  } else return( (edge[1]-1)*n - sum(1:(edge[1]-1)) + (edge[2] - edge[1]) )
}

edge_index_inv = function(index,n) {
  # the first edge is at least index/n
  ind1 = as.integer(max(round(index/n,0),1))
  while ((ind1 * n - sum(1:ind1)) < index) {ind1 = ind1 + 1}
  ind2 = n - ((ind1 * n - sum(1:ind1)) - index) 
  return(as.character(paste(ind1,"-",ind2, sep="")))
}



fast_ewens = function(y, theta_hat) {
  beta = log(theta_hat)
  eta = rep(1,nt) * beta
  mu <- sapply(eta, b1)
  W <- sapply(eta, b2)
  z <- eta + (y - mu) / W 
  wlm <- lm(z ~ 1, weights = W)
  return(coef(wlm))
}
  


# input: g - chordal graph
# output: edges/ added - removed
lazy_proposal = function(g,n) {
  ep = list()
  added = c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if(get.edge.ids(g,c(i,j)) && is_chordal(g-E(g,c(i,j)))$chordal) {
        ep = list.append(ep,c(i,j))
        added = c(added,0)
      } else if( (!(get.edge.ids(g,c(i,j)))) && is_chordal(g %>% add_edges(c(i,j)))$chordal ) {
        ep = list.append(ep,c(i,j))
        added = c(added,1)
      }
    }
  }
  return(list(edge_perturb=ep, added=added))
}

####### Graph Estimate ##########


# input: co-occurence matrix/ item pair
# output: 2 x 2 contingency table
cont_table = function(M, items, nt) {
  both = M[items[1], items[2]] # both items in the transaction
  one_only = M[items[1], items[1]] - both # item 1 but not item 2
  two_only = M[items[2], items[2]] - both # item 1 but not item 2
  neither = nt - both - one_only - two_only
  return(array( c(both, one_only, two_only, neither), dim = c(2,2)) )
}

# input: binary data matrix (D)
# output: graph estimate, G hat, based on fishers exact test with bonferroni correction
g_estimate = function(D, pval = .0001, M, nt) {
  n = ncol(D)
  nt = nrow(D)
  M = t(D) %*% D # co-occurence matrix 
  A = array(0, dim=c(n,n)) # adjacency matrix
  threshold = pval/choose(n,2)
  pvals = array(0,dim=c(choose(n,2), 2))
  k = 1 # counter for data frame
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      pv = fisher.test( cont_table(M, c(i,j), nt) )$p.value
      pvals[k,1] = edge_index(c(i,j),n)
      pvals[k,2] = pv
      k = k + 1
      if(pv <= threshold) {
        A[i,j]=1
        A[j,i]=1
      }
    }
  }
  g2 = graph_from_adjacency_matrix(A, mode="undirected")
  chordalg2 = is_chordal(g2, fillin = T, newgraph = T)
  g3 = chordalg2$newgraph # ghat
  
  C <- lapply(cliques(g3), as.vector) # quick check: `table(sapply(C, length))`
  cn <- max(sapply(C, length)) - 1 # clique number of `g` (minus 1)
  
  # graph estimate induces a design matrix
  X <- matrix(0, nrow = length(C), ncol = cn + n) # design
  for (i in seq_along(C)) {
    kc <- length(C[[i]]) - 1
    if (kc > 0) X[i, kc] <- 1
    X[i, cn + C[[i]]] <- 1 / (kc + 1)
  }
  
  return(list(graph=g3, pvals = pvals, X=X, C = C, cn = cn, fill = chordalg2$fillin))
}

############ GAMMA HAT PROCEDURE #############

# input: graph/ n
# output: gamma hat

gamma_hat = function(graph, n) {
  # initialize design matrix & response vector
  design = array(0, dim = c(choose(n,2),n))
  Y = array(0, dim = c(choose(n,2),1))
  
  # build design matrix & response vector
  row_num = 1
  for (i in 1:(n-1)) {
    for(j in (i+1):n) {
      design[row_num, c(i,j)] = 1
      if (get.edge.ids(graph, c(i,j)) > 0) {Y[row_num] = 1}
      row_num = row_num + 1
    }
  }
  
  # estimate gamma 
  #mod_gamma = glm(Y ~ design-1, family = binomial)
  mod_gamma = bayesglm(Y ~ design-1, family = binomial)
  return(mod_gamma)
}


################################# OPTIMIZATION PROCEDURE FUNCTIONS ##################################

# input: old design/ clique set matrix
# output: new design/ clique set matrix
# newX = function(X, ec, added, ghat, cn=cn_hat) {
#   if(added==0) {
#     del_rows = which(apply(X,1,function(x) (x[ec[1]+cn]>0 && x[ec[2]+cn]>0) ))
#     if (length(del_rows) != 0) X = X[-del_rows,]
#     return(X)
#   } else {
#     grp = c(ec,intersect(neighbors(ghat,ec[1]), neighbors(ghat,ec[2])))
#     gt = induced.subgraph(ghat, grp)
#     V(gt)$name = grp[order(grp)]
#     j = cliques(gt)
#     to_add = j[unlist(lapply(j, function(x) is.subset(ec, as.numeric(labels(x))) ))]
#     
#     Xnew <- matrix(0, nrow = length(to_add), ncol = cn + n) # design
#     for (i in seq_along(to_add)) {
#       kc <- length(to_add[[i]]) - 1
#       if (kc > 0) Xnew[i, kc] <- 1
#       Xnew[i, cn + as.numeric(labels(to_add[[i]])) ] <- 1 /((kc + 1))
#     }
#     return(rbind(X,Xnew))
#   }
# }


getDesign = function(g) {
  C <- lapply(cliques(g), as.vector) 
  cn <- max(sapply(C, length)) - 1
  X <- matrix(0, nrow = length(C), ncol = cn + n) # design
  for (i in seq_along(C)) {
    kc <- length(C[[i]]) - 1
    if (kc > 0) X[i, kc] <- 1
    X[i, cn + C[[i]]] <- 1 /((kc + 1))
  }
  return(X)
}


newX = function(X, g, ec, added, cn=cn_hat) {
    if(added==0) {
      del_rows = which(apply(X,1,function(x) (x[ec[1]+cn]>0 && x[ec[2]+cn]>0) ))
      if (length(del_rows) != 0) X = X[-del_rows,]
      return(X)
    } else {
    gnew = add_edges(g, ec)
    C <- lapply(cliques(gnew), as.vector) 
    cn <- max(sapply(C, length)) - 1
    
    # design matrix for initial graph estimate
    X <- matrix(0, nrow = length(C), ncol = cn + n) # design
    for (i in seq_along(C)) {
      kc <- length(C[[i]]) - 1
      if (kc > 0) X[i, kc] <- 1
      X[i, cn + C[[i]]] <- 1 /((kc + 1))
    }
    return(X)
    }
}

# allow the new S matrix to expand if the clique number of the new graph is bigger than the old graph
# new estimate of S and k given an edge change and chordal graph
sK_update = function(D,S,k,edge_change,g,added,cn) {
  S_new = S
  k_new = k
  # transactions in which the two vertices involved in the edge change co-occur
  impacted = which(apply(D,1,function(x) (x[ edge_change[1] ]>0 && x[ edge_change[2] ] >0 )) > 0)
  
  # if no transactions are impacted, return current S,k
  if (length(impacted) > 0) {
    #print(paste(c(edge_change, length(impacted))))
    if (added==1) {g4=add_edges(g,edge_change)} else {g4 = g - E(g,edge_change)}
    
    # get clique number of new graph
    Cl <- lapply(cliques(g4), as.vector) # quick check: `table(sapply(C, length))`
    cn_g4 <- max(sapply(Cl, length)) - 1 # clique number of `g` (minus 1)
    
    # if the new clique number is greater than the old clique number, the S matrix needs to expand
    # in particular, add 0 columns
    if(cn_g4 > cn) {S_new = cbind(S[,1:cn],matrix(0, nrow = dim(S)[1], ncol = (cn_g4 - cn)),S[,(cn+1):dim(S)[2] ])}
    
    # if the new clique number is LOWER than the old clique number, the S matrix needs to contract
    if(cn_g4 < cn) {S_new = cbind(S[,1:cn_g4],S[,(cn+1):dim(S)[2] ])}
    
    
    kS_new = kestimate(g4,D[impacted,],cn_g4)
    
    k_new[impacted] = kS_new$k
    S_new[impacted,] = kS_new$S
  }
  return(list(k=k_new, S=S_new))
}


# input: chordal g, edge change, current parameter estimates (S,k,alpha, beta, theta)
# output: new parameter estimates (S,k,alpha, beta, theta)
update_params = function(g, ec, added, parameters, M) {
  Xp = newX(parameters$X, g, ec, added) # can make this function more efficient
  
  if(dim(Xp)[2] - n > cn_hat) {Xp = Xp[,-c(cn_hat+1)]}
  
  ng = log(sum(exp(Xp %*% parameters$beta)))
  
  # fast alpha/ beta
  cj = t(apply(Xp,1, function(x) ifelse(x>0,1,0)))
  CO = t(cj) %*% cj
  PI = make_pi(parameters$beta, Xp, cj, CO)
  W = make_W(parameters$beta, PI, Xp, CO) + diag(rep(.001,dim(PI)[1]))
  U = chol(W)
  Uinv = solve(U)
  Winv = Uinv %*% t(Uinv)
  
  if (M[edge_change[1], edge_change[2]] > 1) {
    kS = sK_update(parameters$S, parameters$k, ec, g, added)
    
    # fast theta hat (don't need a new theta hat if there are insufficient co-occurences)
    th = exp(fast_ewens(kS$k, parameters$theta)) 
    
    # new beta in the co-occur case
    beta_tilde = parameters$beta + sum(kS$k)^-1 * Winv %*% (colSums(kS$S)-sum(kS$k)*diag(PI))
    
    # FIX ALPHA
    #beta_tilde[1:cn] = 1:cn
    
    return(list(X = Xp, S=kS$S, k=kS$k, beta=beta_tilde, theta=th, ng=ng))
  } 
  else {
    beta_tilde = parameters$beta + sum(k_hat)^-1 * Winv %*% (colSums(parameters$S)-sum(parameters$k)*diag(PI))
    
    # FIX ALPHA
    #beta_tilde[1:cn] = 1:cn
    
    return(list(X = Xp, S=parameters$S, k=parameters$k, beta=beta_tilde, theta=parameters$theta, ng=ng))
  } 
}






