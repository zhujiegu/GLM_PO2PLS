# this method fit 1 joint-comp glm-PO2PLS binary model on filtered x and y sequentially, 
# updating only one pair of regression coefficients 
Su_PO2PLS_bi_filter <- function(X, Y, Z, r, rx, ry, steps, tol, level, 
                              Nr.core =1, init_param= 'o2m', orth_type = 'SVD', 
                              random_restart = 'F'){
  if(r==1){
    return(Su_PO2PLS_bi(X, Y, Z, 1, rx, ry, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                        init_param= init_param, orth_type = orth_type, random_restart = random_restart))
  }
  if(r>1){
    print('More than 1 joint comp, fitting GLM-PO2PLS Gaussian model first')
    fit_tmp <- Su_PO2PLS(X, Y, Z, r, rx, ry, steps = 1000, tol = tol, init_param= init_param)
    fit_tmp$params <- fit_tmp$params[names(fit_tmp$params)!='sig2G']
    
    for(k in 1:r){
      # filter data
      x_k <- with(fit_tmp, X - latent_var$mu_T[,-k,drop=F] %*% t(params$W[,-k,drop=F]) -
                    latent_var$mu_To %*% t(params$Wo))
      y_k <- with(fit_tmp, Y - latent_var$mu_U[,-k,drop=F] %*% t(params$C[,-k,drop=F]) -
                    latent_var$mu_Uo %*% t(params$Co))
      print(paste('fitting binary model for component',k))
      if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a0","a","b") %in% names(init_param))){
        message('using modified old fit as initial parameters \n')
        params_k <- init_param
        params_k$W <- params_k$W[,k,drop=F]
        params_k$C <- params_k$C[,k,drop=F]
        params_k$B <- params_k$B[k,k,drop=F]
        params_k$SigT <- params_k$SigT[k,k,drop=F] 
        params_k$SigH <- params_k$SigH[k,k,drop=F]
        params_k$SigU <- params_k$SigU[k,k,drop=F]
        params_k$a0 <- params_k$a0
        params_k$a <- params_k$a[,k,drop=F]
        params_k$b <- params_k$b[,k,drop=F]
        params_k$Wo <- params_k$Wo*0
        params_k$Co <- params_k$Co*0
        params_k$SigTo <- matrix(0,1,1)
        params_k$SigUo <- matrix(0,1,1)
        
        fit_k <- Su_PO2PLS_bi(x_k, y_k, Z, 1, 0, 0, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                              init_param= params_k, orth_type = orth_type, random_restart = random_restart)
      }else{
        fit_k <- Su_PO2PLS_bi(x_k, y_k, Z, 1, 0, 0, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                              init_param= init_param, orth_type = orth_type, random_restart = random_restart)
      }
      
      # update initial estimate
      fit_tmp$params$W[,k] <- fit_k$params$W
      fit_tmp$params$C[,k] <- fit_k$params$C
      fit_tmp$params$B[k,k] <- fit_k$params$B %>% as.numeric()
      fit_tmp$params$SigT[k,k] <- fit_k$params$SigT %>% as.numeric()
      fit_tmp$params$SigH[k,k] <- fit_k$params$SigH %>% as.numeric()
      fit_tmp$params$SigU[k,k] <- fit_k$params$SigU %>% as.numeric()
      fit_tmp$params$a0 <- fit_k$params$a0
      fit_tmp$params$a[,k] <- fit_k$params$a
      fit_tmp$params$b[,k] <- fit_k$params$b
      fit_tmp$latent_var$mu_T[,k] <- fit_k$latent_var$mu_T
      fit_tmp$latent_var$mu_U[,k] <- fit_k$latent_var$mu_U
    }
    return(fit_tmp)
  }
}


# A one-comp binary model that updates all the (a,b)
Su_PO2PLS_bi_fast <- function(X, Y, Z, r, rx, ry, steps = 1e5, tol = 1e-6, level=level, 
                              Nr.core =1, init_param= c('random','o2m'), orth_type = "SVD", 
                              random_restart = FALSE){

  # initial fit of GLM-PO2PLS Gaussian modle
  message('Initial fit with GLM-PO2PLS Gaussian model')
  fit_ini <- Su_PO2PLS(X, Y, Z, r, rx, ry, steps = 1000, tol = tol, init_param= init_param)
  fit_ini$params <- fit_ini$params[names(fit_ini$params)!='sig2G']
  
  if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a0","a","b") 
         %in% names(init_param))){
    fit_ini$params$a <- init_param$a
    fit_ini$params$b <- init_param$b
    fit_ini$params$a0 <- init_param$a0
    }else{
      fit_ini$params$a0 <- 0
    }
  params <- fit_ini$params
  fixtu <- with(fit_ini$latent_var, cbind(mu_T, mu_U))
  
  for(k in 1:r){
    # filter data
    x_k <- X - fixtu[, -c(k,r+1:r),drop=F] %*% t(params$W[,-k,drop=F]) -
      fit_ini$latent_var$mu_To %*% t(params$Wo)
    y_k <- Y - fixtu[, -c(k+r,1:r),drop=F] %*% t(params$C[,-k,drop=F]) -
      fit_ini$latent_var$mu_Uo %*% t(params$Co)
    print(paste('fitting binary model for component',k))
    
    logl = 0*0:steps
    tic <- proc.time()
    print(paste('started',date()))
    
    i_rr <- 0
    random_restart_original <- random_restart
    random_restart <- TRUE
    while(random_restart){
      
      if(i_rr > 0) {
        message("Log-likelihood: ", logl[i+1])
        message(paste("random restart no",i_rr))
      }
      params_max <- params
      for(i in 1:steps){
        E_next = E_step_bi_fast(x_k, y_k, Z, params, fixtu, k = k, level=level, Nr.core=Nr.core)
        params_next = M_step_bi(E_next, params, X, Y, Z, orth_type = orth_type)
        params_next$B <- abs(params_next$B)
        
        if(i == 1) logl[1] = E_next$logl
        logl[i+1] = E_next$logl
        
        # # influence of Z
        # z_inf[i] <- with(params_next, sd(X%*%W/sig2E)/sd(Z%*%a/sig2G))
        
        # #debug
        # if(is.na(logl[i+1])) browser()
        # # if(logl[i+1] < logl[i]) browser()
        # if(logl[i+1] >1000) browser()
        
        if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
        if(i %in% c(1, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)) {
          print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
        }
        # solve(0)
        if(logl[i+1] > max(logl[1:i])) params_max <- params_next
        params = params_next
        
        # debug_list[[i]] <- list(Cor=cor(cbind(E_next$mu_T,E_next$mu_U)),
        #                         aa=params_next$a,
        #                         bb=params_next$b,
        #                         sig2G=params_next$sig2G,
        #                         Sttuu=with(E_next, rbind(cbind(2*Stt*diag(1,r), Stu+t(Stu)),
        #                                        cbind(Stu+t(Stu), 2*Suu*diag(1,r)))),
        #                         ZTU=t(Z/N)%*%cbind(E_next$mu_T,E_next$mu_U))
      }
      if(!any(diff(logl[-1]) < -1e-10) | !random_restart_original) {
        random_restart = FALSE
        break
      }
      i_rr <- i_rr + 1
      params <- jitter_params(params)
      params[-(1:4)] <- generate_params_bi(X, Y, Z, r, rx, ry, type = 'r')[-(1:4)]
    }
  }
}


E_step_bi_fast <- function(X, Y, Z, params, fixtu, k, level = level, Nr.core=1){
  # retrieve parameters
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
  # define dimensions
  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(a)
  
  # Numerical integration of the common part and store the results (all sample combined)
  common <- GH_com_fast(Nr.cores=Nr.core, level=level, X=X,Y=Y,Z=Z, params=params, 
                        fixtu=fixtu, k=k, plot_nodes=F)
  list_com_log <- common$list_com_log # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- common$nodes
  # divide both numerator and denominator by exp(max), to avoid "Inf"
  list_com_sub <- lapply(list_com_log, function(e){
    e <- unlist(e)
    max_e <- max(e)
    as.list(exp(e-max_e))
  }
  )
  
  # likelihood of each sample
  list_lik_log <- lapply(list_com_log, function(e){
    e <- unlist(e)
    max_e <- max(e)
    return(max_e + log(sum(exp(e-max_e))))
  })
  
  # Denominator = f(x,y,z)/exp(max)
  dnmt <- lapply(1:N, function(e) Reduce("+", list_com_sub[[e]]))
  
  # E(tu|xyz) with numerical integration
  ## different part in integrand
  int_diff <- list_nodes
  ## Numerator
  mu_TU_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff, list_com_sub[[e]])))
  mu_TU <- Map("/", mu_TU_nu, dnmt)
  mu_TU <- mu_TU %>% unlist %>% matrix(nrow=2) %>% t
  
  # E(Sttuu|xyz) with numerical integration
  ## different part in integrand
  int_diff <- lapply(list_nodes, crossprod)
  ## Numerator
  S_ttuu_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff, list_com_sub[[e]])))
  S_ttuu <- Map("/", S_ttuu_nu, dnmt)
  # 2r x 2r, Average across N samples
  # S_ttuu <- Reduce("+", S_ttuu)/N
  # Stt <- S_ttuu[1:r,1:r,drop=FALSE]
  # Suu <- S_ttuu[r+1:r,r+1:r,drop=FALSE]
  # Stu <- S_ttuu[1:r,r+1:r,drop=FALSE]
  # Chh <- Suu - t(Stu)%*%B - t(B)%*%Stu + t(B)%*%Stt%*%B
  
  # adjust (t,u) to (t,h)
  list_nodes_th <- lapply(list_nodes, function(e) {
    e[,r+1:r] <- e[,r+1:r] - e[,1:r]%*% B
    return(e)
  })
  beta <- cbind(a0,a,b)
  
  # Numerical estimation of Q_ab and gradient
  Q_ab <- GH_Q_ab(beta,list_nodes_th,Z,list_com_sub,dnmt)
  grd_ab <- GH_grd_ab(beta,list_nodes_th,Z,list_com_sub,dnmt)
  # Backtracking rule to find the step size s
  s = 1
  Q_ab_new <- GH_Q_ab(beta + s*grd_ab,list_nodes_th,Z,list_com_sub,dnmt)
  
  while(Q_ab_new < Q_ab + 0.5*s*tcrossprod(grd_ab)){
    s = 0.8*s
    # print(s)
    Q_ab_new <- GH_Q_ab(beta + s*grd_ab,list_nodes_th,Z,list_com_sub,dnmt)
  }
  
  # log likelihood
  logl <- list_lik_log %>% unlist %>% sum
  # print(list_lik %>% unlist %>% log)
  
  # output
  list(
    mu_T = mu_TU[,1:r,drop=F],
    mu_U = mu_TU[,r+1:r,drop=F],
    S_ttuu = S_ttuu,
    # Stt = Stt,
    # Suu = Suu,
    # Stu = Stu,
    # Shh = Chh,
    s = s,
    grd_ab = grd_ab,
    logl = logl,
    GH_common=common
  )
}


# common parts in GH for one component model
GH_com_fast <- function(Nr.cores=1, level=5, X,Y,Z, params, fixtu, k, plot_nodes=F){
  dim = 2
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
  sigma <- with(params, rbind(cbind(SigT[k,k,drop=F], SigT[k,k,drop=F]%*%B[k,k,drop=F]),
                              cbind(B[k,k,drop=F]%*%SigT[k,k,drop=F], 
                                    SigT[k,k,drop=F]%*%B[k,k,drop=F]^2 + SigH[k,k,drop=F])))
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
  
  list_com_log <- parallel::mclapply(1:N, mc.cores = Nr.cores, function(e){
    lapply(n_w, fun_com_fast, x=X[e,],y=Y[e,],z=Z[e], params=params, k=k, fixtu_i=fixtu[e,,drop=F])
  })
  # browser()
  
  # lapply(n_w, fun_com, x=x[e,],y=y[e,],z=z[e], params=params)
  
  nodes <- lapply(seq_len(nrow(nodes)), function(i) nodes[i, ,drop=F])
  return(list(list_com_log=list_com_log, nodes=nodes))
}


# common part of every integral
# f(z|tu) * f(x|t) * f(y|u)
fun_com_fast <- function(i, x,y,z,params,k,fixtu_i){
  r <- 1
  rr <- ncol(fixtu)/2
  p <- nrow(params$W)
  q <- nrow(params$C)
  alpha <- with(params, cbind(a0,a-b%*%B,b))
  fixtu_i[ ,c(k, rr+k)] <- i$n
  z_tu <- with(params, 
               if(z==1){
                 1/(1+exp(-(cbind(1,fixtu_i) %*% t(alpha))))
               }else{
                 1/(1+exp((cbind(1,fixtu_i) %*% t(alpha))))
               }
  ) %>% as.numeric()

  # log(f(x|t))
  l_xt <- with(params, -0.5*(
    p*log(2*pi*sig2E) + 1/sig2E*(
      (x-i$n[,1]%*%t(W[,k]))%*%t(x-i$n[,1]%*%t(W[,k]))))) %>% as.numeric()
  
  # log(f(y|u))
  l_yu <- with(params, -0.5*(
    q*log(2*pi*sig2F) + 1/sig2F*(
      (y-i$n[,2]%*%t(C[,k]))%*%t(y-i$n[,2]%*%t(C[,k]))))) %>% as.numeric()
  
  return(log(z_tu) + l_xt + l_yu + log(i$w))
}
