GH_Intl <- function(func, dim=2*r, level=6, x,y,z, params, plot_nodes=F){
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

  # calculate each node
  tmp <- lapply(n_w, func, x=x,y=y,z=z, params=params)
  # print(tmp[[1]])
  # tmp <- lapply(1:level^dim, fun_h, nodes=nodes, w=w)
  # sum over level^dim combos
  result <- Reduce("+", tmp)
  # print(result)
  
  # need fixing here
  
  tmp_com <- lapply(n_w, fun_com, x=x,y=y,z=z, params=params)
  # tmp <- lapply(1:level^dim, fun_h, nodes=nodes, w=w)
  # sum over level^dim combos
  result_com <- Reduce("+", tmp_com)
  # print(result_com)
  
  return(result/result_com)
}

# # checking integral with covariance matrix
# fun_h <- function(i, x,y,z,params){
#   crossprod(i$n) * i$w
# }

# E(tu|xyz)
fun_mu <- function(i, x,y,z,params){
  tu <- i$n
  z_tu <- with(params, 
               if(z==1){
                 1/(1+exp(-(cbind(1,i$n) %*% t(params$alpha))))
               }else{
                 1/(1+exp((cbind(1,i$n) %*% t(params$alpha))))
               }
  ) %>% as.numeric()
  
  Sig_xt <- with(params, Wo %*% SigTo %*% t(Wo) + diag(sig2E,p))
  x_t <- with(params,
              (2*pi)^(-0.5*p) * det(Sig_xt)^(-0.5) *
                exp(-0.5*(x-i$n[,1:r]%*%t(W))%*%solve(Sig_xt)%*%t((x-i$n[,1:r]%*%t(W))))
  ) %>% as.numeric()
  
  Sig_yu <- with(params, Co %*% SigUo %*% t(Co) + diag(sig2F,q))
  y_u <- with(params,
              (2*pi)^(-0.5*q) * det(Sig_yu)^(-0.5) *
                exp(-0.5*(y-i$n[,r+1:r]%*%t(C))%*%solve(Sig_yu)%*%t((y-i$n[,r+1:r]%*%t(C))))
  ) %>% as.numeric()
  return(tu *z_tu*x_t*y_u* i$w)
}

# E(Sttuu|xyz)
fun_S <- function(i, x,y,z,params){
  Sttuu <- crossprod(i$n)
  z_tu <- with(params, 
               if(z==1){
                 1/(1+exp(-(cbind(1,i$n) %*% t(params$alpha))))
               }else{
                 1/(1+exp((cbind(1,i$n) %*% t(params$alpha))))
               }
  ) %>% as.numeric()
  
  Sig_xt <- with(params, Wo %*% SigTo %*% t(Wo) + diag(sig2E,p))
  x_t <- with(params,
              (2*pi)^(-0.5*p) * det(Sig_xt)^(-0.5) *
                exp(-0.5*(x-i$n[,1:r]%*%t(W))%*%solve(Sig_xt)%*%t((x-i$n[,1:r]%*%t(W))))
  ) %>% as.numeric()
  
  Sig_yu <- with(params, Co %*% SigUo %*% t(Co) + diag(sig2F,q))
  y_u <- with(params,
              (2*pi)^(-0.5*q) * det(Sig_yu)^(-0.5) *
                exp(-0.5*(y-i$n[,r+1:r]%*%t(C))%*%solve(Sig_yu)%*%t((y-i$n[,r+1:r]%*%t(C))))
  ) %>% as.numeric()
  return(Sttuu*z_tu*x_t*y_u*i$w)
}


#################
# isolates this part for computational efficiency later
# f(z|tu) * f(x|t) * f(y|u)
fun_com <- function(i, x,y,z,params){
  z_tu <- with(params, 
               if(z==1){
                 1/(1+exp(-(cbind(1,i$n) %*% t(params$alpha))))
               }else{
                 1/(1+exp((cbind(1,i$n) %*% t(params$alpha))))
               }
               ) %>% as.numeric()
  
  Sig_xt <- with(params, Wo %*% SigTo %*% t(Wo) + diag(sig2E,p))
  x_t <- with(params,
              (2*pi)^(-0.5*p) * det(Sig_xt)^(-0.5) *
                exp(-0.5*(x-i$n[,1:r]%*%t(W))%*%solve(Sig_xt)%*%t((x-i$n[,1:r]%*%t(W))))
  ) %>% as.numeric()
  
  Sig_yu <- with(params, Co %*% SigUo %*% t(Co) + diag(sig2F,q))
  y_u <- with(params,
              (2*pi)^(-0.5*q) * det(Sig_yu)^(-0.5) *
                exp(-0.5*(y-i$n[,r+1:r]%*%t(C))%*%solve(Sig_yu)%*%t((y-i$n[,r+1:r]%*%t(C))))
  ) %>% as.numeric()
  return(z_tu*x_t*y_u*i$w)
}