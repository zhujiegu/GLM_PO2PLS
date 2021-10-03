#' @export
E_step_bi <- function(X, Y, Z, params, level = level, Nr.core=1){
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
  
  #################################################
  ## Numerical integration of the common part and store the results (all sample combined)
  common <- GH_com(Nr.cores=Nr.core, level=level, X,Y,Z, params, plot_nodes=F)
  list_com_log <- common$list_com_log # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- common$nodes
  # divide both numerator and denominator by exp(max), to avoid "Inf"
  list_com_sub <- lapply(list_com_log, function(e){
    e <- unlist(e)
    max_e <- max(e)
    as.list(exp(e-max_e))
  }
  )
  #################################################
  ## likelihood of each sample
  list_lik_log <- lapply(list_com_log, function(e){
    e <- unlist(e)
    max_e <- max(e)
    return(max_e + log(sum(exp(e-max_e))))
  })
  #################################################
  ## when likelihood is numerically 0
  # idx_zero <- which(unlist(list_lik) == 0)
  #################################################
  # Denominator = f(x,y,z)/exp(max)
  dnmt <- lapply(1:N, function(e) Reduce("+", list_com_sub[[e]]))
  #################################################
  ## # E(tu|xyz) with numerical integration
  # different part in integrand
  int_diff <- list_nodes
  # Numerator
  mu_TU_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff, list_com_sub[[e]])))
  mu_TU <- Map("/", mu_TU_nu, dnmt)
  mu_TU <- mu_TU %>% unlist %>% matrix(nrow=2) %>% t
  #################################################
  ## # E(Sttuu|xyz) with numerical integration
  # different part in integrand
  int_diff <- lapply(list_nodes, crossprod)
  # Numerator
  S_ttuu_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff, list_com_sub[[e]])))
  S_ttuu <- Map("/", S_ttuu_nu, dnmt)
  
  # # 2r x 2r, Average across N samples
  # S_ttuu <- Reduce("+", S_ttuu)/N
  
  # Stt <- S_ttuu[1:r,1:r,drop=FALSE]
  # Suu <- S_ttuu[r+1:r,r+1:r,drop=FALSE]
  # Stu <- S_ttuu[1:r,r+1:r,drop=FALSE]
  # Chh <- Suu - t(Stu)%*%B - t(B)%*%Stu + t(B)%*%Stt%*%B
  #################################################
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
  
  #################################################
  ## log likelihood
  logl <- list_lik_log %>% unlist %>% sum
  # print(list_lik %>% unlist %>% log)
  #################################################
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
    logl = logl
  )
}

## function for numerical integration of Q_ab
GH_Q_ab <- function(beta, l_n=list_nodes_th,Z.=Z,l_com=list_com_sub,de=dnmt){
  # log probability of z=1 and z=0
  z1 <- lapply(l_n, function(e) log(as.numeric(1/(1+exp(-(cbind(1,e) %*% t(beta)))))))
  z0 <- lapply(l_n, function(e) log(as.numeric(1/(1+exp((cbind(1,e) %*% t(beta)))))))
  # different part in integrand
  int_diff_z <- lapply(Z., function(z){
    if(z){
      return(z1)
    }else{
      return(z0)
    }
  })
  
  # combine common parts
  Q_ab <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff_z[[e]], l_com[[e]])))
  Q_ab <- Map("/", Q_ab, de)
  # Add N samples
  Q_ab <- Reduce("+", Q_ab)
  return(Q_ab)
}

## function for numerical integration of gradient of Q_ab
GH_grd_ab <- function(beta, l_n=list_nodes_th,Z.=Z,l_com=list_com_sub,de=dnmt){
  # Q'_ab
  z1 <- lapply(l_n, function(e) {
    cbind(1,e)*as.numeric(exp(-(cbind(1,e) %*% t(beta))))/
      (1+as.numeric(exp(-(cbind(1,e) %*% t(beta)))))
  })
  z0 <- lapply(l_n, function(e) {
    -cbind(1,e)*as.numeric(exp((cbind(1,e) %*% t(beta))))/
      (1+as.numeric(exp((cbind(1,e) %*% t(beta)))))
  })
  # different part in integrand
  int_diff_z <- lapply(Z., function(z){
    if(z){
      return(z1)
    }else{
      return(z0)
    }
  })
  # combine common parts
  grd_ab <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff_z[[e]], l_com[[e]])))
  grd_ab <- Map("/", grd_ab, de)
  # Add N samples
  grd_ab <- Reduce("+", grd_ab)
  return(grd_ab)
}

#' @export
M_step_bi <- function(E_fit, params, X, Y, Z, orth_type = c("SVD","QR")){
  orth_x = ssq(params$Wo) > 0
  orth_y = ssq(params$Co) > 0
  orth_type = match.arg(orth_type)
  with(E_fit,{
    
    N = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    r = ncol(mu_T)
    rx = ncol(params$Wo)
    ry = ncol(params$Co)
    
    # 2r x 2r, Average across N samples
    S_ttuu_avgN <- Reduce("+", S_ttuu)/N
    Stt <- S_ttuu_avgN[1:r,1:r,drop=FALSE]
    Suu <- S_ttuu_avgN[r+1:r,r+1:r,drop=FALSE]
    Stu <- S_ttuu_avgN[1:r,r+1:r,drop=FALSE]
    Shh <- Suu - t(Stu)%*%params$B - t(params$B)%*%Stu + t(params$B)%*%Stt%*%params$B
    
    # print(S_ttuu_avgN)
    # print(cbind(Stt, Suu, Shh))
    
    
    # filter out samples with negative sig2E or sig2F
    # Stt, Suu different for each sample
    tmp_sig2E = sapply(1:N, function(e){
      1/p * (sum(diag(X[e,,drop=F]%*%t(X[e,,drop=F]))) - 2*sum(diag(mu_T[e,,drop=F]%*%t(params$W)%*%t(X[e,,drop=F]))) +
               sum(diag(S_ttuu[[e]][1:r,1:r,drop=FALSE])) - sum(diag(params$SigTo)))
    })
    tmp_sig2F = sapply(1:N, function(e){
      1/q * (sum(diag(Y[e,,drop=F]%*%t(Y[e,,drop=F]))) - 2*sum(diag(mu_U[e,,drop=F]%*%t(params$C)%*%t(Y[e,,drop=F]))) +
               sum(diag(S_ttuu[[e]][r+1:r,r+1:r,drop=FALSE])) - sum(diag(params$SigUo)))
    })
    
    # # First average Stt, Suu, same for all samples
    # tmp_sig2E = sapply(1:N, function(e){
    #   1/p * (sum(diag(X[e,,drop=F]%*%t(X[e,,drop=F]))) - 2*sum(diag(mu_T[e,,drop=F]%*%t(params$W)%*%t(X[e,,drop=F]))) +
    #            sum(diag(Stt)) - sum(diag(params$SigTo)))
    # })
    # tmp_sig2F = sapply(1:N, function(e){
    #   1/q * (sum(diag(Y[e,,drop=F]%*%t(Y[e,,drop=F]))) - 2*sum(diag(mu_U[e,,drop=F]%*%t(params$C)%*%t(Y[e,,drop=F]))) +
    #            sum(diag(Suu)) - sum(diag(params$SigUo)))
    # })
    
    # # which samples have positive sig2E or sig2F
    # idx_posE <- which(tmp_sig2E>0)
    # idx_posF <- which(tmp_sig2F>0)
    
    # if(length(idx_posE)==0 | length(idx_posF)==0){
    #   stop("all samples have negative sig2E or sig2F")
    # }
    
    # nr_negE <- N - length(idx_posE)
    # nr_negF <- N - length(idx_posF)
    
    # print(nr_negE); print(nr_negF)
    # print(tmp_sig2E)
    # print(tmp_sig2F)
    
    params$sig2E = mean(tmp_sig2E)
    params$sig2F = mean(tmp_sig2F)
    
    params$B = t(Stu) %*% MASS::ginv(Stt) * diag(1,r)
    params$SigT = Stt*diag(1,r)
    params$SigU = Suu*diag(1,r)
    params$SigH = Shh*diag(1,r)
    params$a0 = params$a0 + s*grd_ab[,1]
    params$a = params$a + s*grd_ab[,1+1:r]
    params$b = params$b + s*grd_ab[,1+r+1:r]
    
    params$W = orth(t(X/N) %*% mu_T %*% MASS::ginv(Stt),type = orth_type)
    params$C = orth(t(Y/N) %*% mu_U %*% MASS::ginv(Suu),type = orth_type)#
    
    # params$sig2E = 1/p * abs(sum(diag(X%*%t(X)))/N - 2*sum(diag(mu_T%*%t(params$W)%*%t(X)))/N +
    #                        sum(diag(params$SigT)) - sum(diag(params$SigTo)))
    # params$sig2F = 1/q * abs(sum(diag(Y%*%t(Y)))/N - 2*sum(diag(mu_U%*%t(params$C)%*%t(Y)))/N +
    #                        sum(diag(params$SigU)) - sum(diag(params$SigUo)))
    
    # SVD for high dimensional
    dcmp_varxo <- svd_orthpart(X, params$W, params$sig2E, mu_T, Stt, rx)
    dcmp_varyo <- svd_orthpart(Y, params$C, params$sig2F, mu_U, Suu, ry)

    params$Wo = dcmp_varxo$V
    params$SigTo = diag(x=dcmp_varxo$D, nrow=rx)
    params$Co = dcmp_varyo$V
    params$SigUo = diag(x=dcmp_varyo$D, nrow=ry)

    # ## this only works for low dimension, use power iteration for high dimension
    # var_xo <- t(X)%*%X/N -  t(X)%*%mu_T%*%t(params$W)/N - t(t(X)%*%mu_T%*%t(params$W))/N +
    #   params$W%*%Stt%*%t(params$W) - diag(params$sig2E, p)
    # var_yo <- t(Y)%*%Y/N -  t(Y)%*%mu_U%*%t(params$C)/N - t(t(Y)%*%mu_U%*%t(params$C))/N +
    #   params$C%*%Suu%*%t(params$C) - diag(params$sig2F, q)
    # 
    # dcmp_varxo <- svd(var_xo, nu = min(N, rx), nv=0)
    # dcmp_varyo <- svd(var_yo, nu = min(N, ry), nv=0)
    # params$Wo = dcmp_varxo$u
    # params$SigTo = diag(x=dcmp_varxo$d[1:rx], nrow=rx)
    # params$Co = dcmp_varyo$u
    # params$SigUo = diag(x=dcmp_varyo$d[1:ry], nrow=ry)
    
    # print(str(params))

    
    # Orthogonal loadings are eigen vectors, variance matrices contain the eigen values

    return(params)
  })
}

svd_orthpart <- function(dat, ld, sig, mu, S, rs, tol=1e-10, max_powerit = 10000){
  N <- nrow(dat)
  p <- ncol(dat)
  V <- matrix(NA, nrow = p, ncol = rs)
  D <- c()
  for(comp in 1:rs){
    if(comp == 1){
      # initialize
      v_old <- matrix(rep(1,p),ncol = 1)
      for(i in 1:max_powerit){
        v <- t(dat)%*%(dat%*%v_old)/N - t(dat)%*%mu%*%(t(ld)%*%v_old)/N -
                    ld%*%t(mu)%*%(dat%*%v_old)/N + ld%*%S%*%(t(ld)%*%v_old) - sig*v_old
        d <- norm_vec(v)
        v <- v/d
        if (mse(v_old, v) < tol){
          break
        }
        v_old <- v
      }
      V[,comp] <- v
      D[comp] <- d
    }
    # for high dimension, avoid subtracting pxp matrix 
    if(comp > 1){
      # initialize
      v_old <- matrix(rep(1,p),ncol = 1)
      for(i in 1:max_powerit){
        v <- t(dat)%*%(dat%*%v_old)/N - t(dat)%*%mu%*%(t(ld)%*%v_old)/N -
                    ld%*%t(mu)%*%(dat%*%v_old)/N + ld%*%S%*%(t(ld)%*%v_old) - sig*v_old -
                    V[,1:comp-1]%*%diag(D, nrow = comp-1)%*%(t(V[,1:comp-1])%*%v_old)
        d <- norm_vec(v)
        v <- v/d
        if (mse(v_old, v) < tol){
          break
        }
        v_old <- v
      }
      V[,comp] <- v
      D[comp] <- d
    }
  }
  return(list(V = V, D = D))
}


#' @export
jitter_params <- function(params, amount = NULL){
  suppressWarnings(params[1:4] <- lapply(params[1:4], function(e) sign(ssq(e))*orth(jitter(e,amount = 1))))
  params
}


#' @export
Su_PO2PLS_bi <- function(X, Y, Z, r, rx, ry, steps = 1e5, tol = 1e-6, level=level, Nr.core =1, init_param= c('random','o2m'),
                      orth_type = "SVD", random_restart = FALSE){
  # if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F") %in% names(init_param))) {cat('using old fit \n'); params <- init_param}
  # else {params <- generate_params(X, Y, Z, r, rx, ry, type = init_param)}
  
  # debug_list <- list()
  # z_inf <- c()
  if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a0","a","b") %in% names(init_param))){
    message('using old fit \n')
    params <- init_param}
  else {
    init_param <- match.arg(init_param)
    if(r+max(rx,ry) <= min(ncol(X),ncol(Y)) && init_param == "o2m")
    {
      params <- generate_params_bi(X, Y, Z, r, rx, ry, type = "o2m")
    }
    else
    {
      if(r+max(rx,ry) > min(ncol(X),ncol(Y)) && init_param == "o2m")
      {
        cat("** NOTE: Too many components for init_param='o2m', switched to init_param='unit'**.\n\n");
        init_param = "unit"
      }
      params <- generate_params_bi(X, Y, Z, r, rx, ry, type = init_param)
    }
  }
  
  # params <- generate_params(X, Y, Z, r, rx, ry, type = init_param)
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
      E_next = E_step_bi(X, Y, Z, params, level=level, Nr.core=Nr.core)
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
    params[-(1:4)] <- generate_params_su(X, Y, r, rx, ry, type = 'r')[-(1:4)]
  }
  # params <- params_max
  # signB <- sign(diag(params$B))
  # params$B <- params$B %*% diag(signB,r)
  # params$C <- params$C %*% diag(signB,r)
  # params$b <- params$b %*% diag(signB,r)
  
  ## order of comps for r>1
  # ordSB <- order(diag(params$SigT), decreasing = TRUE)
  # params$B <- params$B[ordSB,ordSB, drop=FALSE]
  # params$W <- params$W[,ordSB, drop=FALSE]
  # params$C <- params$C[,ordSB, drop=FALSE]
  # params$a <- params$a[,ordSB, drop=FALSE]
  # params$b <- params$b[,ordSB, drop=FALSE]
  # params$SigT <- params$SigT[ordSB,ordSB, drop=FALSE]
  # params$SigU <- params$SigU[ordSB,ordSB, drop=FALSE]
  # params$SigH <- params$SigH[ordSB,ordSB, drop=FALSE]
  
  # # R-square Z fit 
  # pre_z <- with(params, X%*%W%*%t(a) + Y%*%C%*%t(b))
  # zfit <- cor(pre_z,Z)^2
  
  message("Nr steps was ", i)
  message("Negative increments: ", any(diff(logl[1:i+1]) < 0), "\n",
          "Smallest increments: ", min(diff(logl[1:i+1])),
          "; Last increment: ", signif(logl[i+1]-logl[i],4))
  message("Log-likelihood: ", logl[i+1])
  # message("R-square of z fit: ", zfit)
  outputt <- list(params = params, logl = logl[0:i+1][-1]) #, debug = debug_list)
  class(outputt) <- "PO2PLS"
  return(outputt)
}