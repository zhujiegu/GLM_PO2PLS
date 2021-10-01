GH_com <- function(Nr.cores=1, level=6, X,Y,Z, params, plot_nodes=F){
  dim = 2*ncol(params$B)
  N = nrow(X)
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
  
  # make list
  n_w <- lapply(seq_len(nrow(nodes)), function(i) list(n = nodes[i, ,drop=F], w = w[i]))
  
  list_com <- mclapply(1:N, mc.cores = Nr.cores, function(e){
    lapply(n_w, fun_com, x=X[e,],y=Y[e,],z=Z[e], params=params)
  })
  # browser()
  
  # lapply(n_w, fun_com, x=x[e,],y=y[e,],z=z[e], params=params)

  nodes <- lapply(seq_len(nrow(nodes)), function(i) nodes[i, ,drop=F])
  return(list(list_com=list_com, nodes=nodes))
}




# func is the g(eta) in form int(g(eta)f(z|eta)f(x,y|eta)f(eta)d(eta))
# set div_mrg to TRUE when divided by f(x,y,z)
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


#################
# common part of every integral
# f(z|tu) * f(x|t) * f(y|u)
fun_com <- function(i, x,y,z,params){
  alpha <- with(params, cbind(a0,a-b%*%B,b))
  z_tu <- with(params, 
               if(z==1){
                 1/(1+exp(-(cbind(1,i$n) %*% t(alpha))))
               }else{
                 1/(1+exp((cbind(1,i$n) %*% t(alpha))))
               }
               ) %>% as.numeric()
  
  # Sig_xt <- with(params, Wo %*% SigTo %*% t(Wo) + diag(sig2E,p))
  rx <- ncol(params$Wo)
  # Sig_xt_inv <- with(params, 1/sig2E*(diag(p) - 1/sig2E*Wo%*%MASS::ginv(MASS::ginv(SigTo)+1/sig2E*diag(rx))%*%t(Wo)))
  # Sig_xt_det <- with(params, det(diag(rx)+1/sig2E*SigTo)*(sig2E^p))
  # x_t_old <- with(params,
  #             (2*pi)^(-0.5*p) * Sig_xt_det^(-0.5) *
  #               exp(-0.5*(x-i$n[,1:r]%*%t(W))%*%MASS::ginv(Sig_xt)%*%t((x-i$n[,1:r]%*%t(W))))
  # ) %>% as.numeric()
  
  # x_t_old <- with(params,
  #             (2*pi)^(-0.5*p) * Sig_xt_det^(-0.5) *
  #               exp(-0.5/sig2E*(
  #                 (x-i$n[,1:r]%*%t(W))%*%t(x-i$n[,1:r]%*%t(W)) -
  #                   1/sig2E*(
  #                     (x-i$n[,1:r]%*%t(W))%*%Wo%*%MASS::ginv(MASS::ginv(SigTo)+1/sig2E*diag(rx))%*%
  #                       (t(Wo)%*%t(x-i$n[,1:r]%*%t(W)))
  #                   )
  #                 ))) %>% as.numeric()
  # print(all.equal(x_t_old,x_t))
  Sig_xt_det_log <- with(params, log(det(diag(rx)+1/sig2E*SigTo)) + p*log(sig2E))
  x_t <- with(params, exp(-0.5*(
    p*log(2*pi) + Sig_xt_det_log + 1/sig2E*(
      (x-i$n[,1:r]%*%t(W))%*%t(x-i$n[,1:r]%*%t(W)) -
        1/sig2E*(
          (x-i$n[,1:r]%*%t(W))%*%Wo%*%MASS::ginv(MASS::ginv(SigTo)+1/sig2E*diag(rx))%*%
            (t(Wo)%*%t(x-i$n[,1:r]%*%t(W)))
        ))))) %>% as.numeric()

  # # without orthogonal part
  # x_t_noot <- with(params, exp(-0.5*(
  #   p*log(2*pi) + p*log(sig2E) +
  #     1/sig2E*((x-i$n[,1:r]%*%t(W))%*%t(x-i$n[,1:r]%*%t(W)))))) %>% as.numeric()
  
  # if(!all.equal(x_t,x_t_old)) stop("x_t wrong")
  # browser()
  if(x_t==Inf) browser()
  
  ry <- ncol(params$Co)
  # Sig_yu <- with(params, Co %*% SigUo %*% t(Co) + diag(sig2F,q))
  # Sig_yu_inv <- with(params, 1/sig2F*(diag(q) - 1/sig2F*Co%*%MASS::ginv(MASS::ginv(SigUo)+1/sig2F*diag(ry))%*%t(Co)))
  # Sig_yu_det <- with(params, det(diag(ry)+1/sig2F*SigUo)*(sig2F^q))
  # y_u <- with(params,
  #             (2*pi)^(-0.5*q) * Sig_yu_det^(-0.5) *
  #               exp(-0.5*(y-i$n[,r+1:r]%*%t(C))%*%MASS::ginv(Sig_yu)%*%t((y-i$n[,r+1:r]%*%t(C))))
  # ) %>% as.numeric()
  
  # y_u <- with(params,
  #             (2*pi)^(-0.5*q) * Sig_yu_det^(-0.5) *
  #               exp(-0.5/sig2F*(
  #                 (y-i$n[,r+1:r]%*%t(C))%*%t(y-i$n[,r+1:r]%*%t(C)) -
  #                   1/sig2F*(
  #                     (y-i$n[,r+1:r]%*%t(C))%*%Co%*%MASS::ginv(MASS::ginv(SigUo)+1/sig2F*diag(ry))%*%
  #                       (t(Co)%*%t(y-i$n[,r+1:r]%*%t(C)))
  #                   )
  #                 ))) %>% as.numeric()
  
  Sig_yu_det_log <- with(params, log(det(diag(ry)+1/sig2F*SigUo)) + q*log(sig2F))
  y_u <- with(params, exp(-0.5*(
    q*log(2*pi) + Sig_yu_det_log + 1/sig2F*(
      (y-i$n[,r+1:r]%*%t(C))%*%t(y-i$n[,r+1:r]%*%t(C)) -
        1/sig2F*(
          (y-i$n[,r+1:r]%*%t(C))%*%Co%*%MASS::ginv(MASS::ginv(SigUo)+1/sig2F*diag(ry))%*%
            (t(Co)%*%t(y-i$n[,r+1:r]%*%t(C)))
        ))))) %>% as.numeric()

  # print(all.equal(y_u_old,y_u))
  # browser()
  # print(c(z_tu, x_t, y_u, i$w))
  ## test likelihood without z
  # return(x_t*y_u*i$w)
  
  return(z_tu*x_t*y_u*i$w)
}