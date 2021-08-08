E_step_logl_bi <- function(X, Y, Z, params, level = level, Nr.core=1){
  ## retrieve parameters
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  SigT = params$SigT
  SigTo = (ssq(Wo)>0)*params$SigTo +0# + 1e-10*(ssq(Wo)==0)
  SigUo = (ssq(Co)>0)*params$SigUo +0# + 1e-10*(ssq(Co)==0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT%*%B^2 + SigH
  a0 = params$a0
  a = params$a
  b = params$b
  ## define dimensions
  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  
  
  # loglikelihood
  loglik <- mclapply(1:N, mc.cores = Nr.core, function(e){
    GH_Intl(fun_1, div_mrg = F, dim=2*r, level=level, X[e,],Y[e,],Z[e], params)}) %>% unlist %>% log %>% sum
  
  return(loglik)
}



GH_Intl <- function(func, div_mrg=F, dim=2*r, level=6, x,y,z, params, plot_nodes=F){
  # standard GH rule
  rule <- fastGHQuad::gaussHermiteData(level)
  
  # expand grid
  nodes <- as.matrix(expand.grid(lapply(apply(replicate(dim, rule$x), 
                                              2, list), unlist)))
  g <- as.matrix(expand.grid(lapply(apply(replicate(dim, rule$w), 
                                          2, list), unlist)))
  w <- apply(g, 1, prod)
  
  # adjust for mu and sigma
  mu <- rep(0,dim)
  sigma <- with(params, rbind(cbind(SigT, SigT%*%B),
                              cbind(B%*%SigT, SigT%*%B^2 + SigH)))
  nodes <- mu + t(sqrt(2)*t(chol(sigma))%*%t(nodes))
  w <- (1/sqrt(pi))^dim * w
  
  # visulize the nodes
  if(plot_nodes){
    plot(nodes, cex=-5/log(w), pch=19,
         xlab=expression(x[1]),
         ylab=expression(x[2]))
  }
  ## plot check with mvQuad package
  # myRule_GH <- function(l){
  #   rule <- fastGHQuad::gaussHermiteData(level)
  #   n <- rule$x
  #   w <- rule$w
  #   initial.domain <- matrix(c(-Inf, Inf), ncol=2)
  #   return(list(n=as.matrix(n), w=as.matrix(w), features=list(initial.domain=initial.domain)))
  # }
  # print(1)
  # nw_myGH <- mvQuad::createNIGrid(d=2, type = "myRule_GH", level = 20)
  # mvQuad::rescale(nw_myGH, m = rep(0,2*r), C = sigma, dec.type = 2)
  # plot(nw_myGH)
  
  # make list
  n_w <- lapply(seq_len(nrow(nodes)), function(i) list(n = nodes[i, ,drop=F], w = w[i]))
  
  # the different part in the integral (g(eta))
  l_part1 <- lapply(n_w, func, x=x,y=y,z=z, params=params)
  
  # the common part in the integral
  l_part2 <- lapply(n_w, fun_com, x=x,y=y,z=z, params=params)
  
  # combine the two parts
  l_result <- Map("*", l_part1, l_part2)
  
  if(div_mrg){
    mrg <- Reduce("+", l_part2)
    result <- Reduce("+", l_result)/mrg
  }else{
    result <- Reduce("+", l_result)
  }
  return(result)
}