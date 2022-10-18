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

# this only update k-th component on the filtered dataset
Su_PO2PLS_bi_filter_k <- function(X, Y, Z, r, rx, ry, k=1, steps, tol, level, 
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
    
    # filter data
    x_k <- with(fit_tmp, X - latent_var$mu_T[,-k,drop=F] %*% t(params$W[,-k,drop=F]) -
                  latent_var$mu_To %*% t(params$Wo))
    y_k <- with(fit_tmp, Y - latent_var$mu_U[,-k,drop=F] %*% t(params$C[,-k,drop=F]) -
                  latent_var$mu_Uo %*% t(params$Co))
    print(paste('fitting binary model for component',k))
    if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a0","a","b")
           %in% names(init_param))){
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
    }else{
      params_k <- fit_tmp$params
      params_k$W <- params_k$W[,k,drop=F]
      params_k$C <- params_k$C[,k,drop=F]
      params_k$B <- params_k$B[k,k,drop=F]
      params_k$SigT <- params_k$SigT[k,k,drop=F] 
      params_k$SigH <- params_k$SigH[k,k,drop=F]
      params_k$SigU <- params_k$SigU[k,k,drop=F]
      params_k$a0 <- 0
      params_k$a <- params_k$a[,k,drop=F]
      params_k$b <- params_k$b[,k,drop=F]
      params_k$Wo <- params_k$Wo*0
      params_k$Co <- params_k$Co*0
      params_k$SigTo <- matrix(0,1,1)
      params_k$SigUo <- matrix(0,1,1)
    }
    fit_k <- Su_PO2PLS_bi(x_k, y_k, Z, 1, 0, 0, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                            init_param= params_k, orth_type = orth_type, random_restart = random_restart)
    return(fit_k)
  }
}

# A one-comp binary model that updates all the (a,b) sequentially for k=1,2,...,K
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
        E_next = E_step_bi_fast( X=x_k, Y=y_k, Z=Z, params=params, fixtu=fixtu, k=k, level=level, 
                                 Nr.core=Nr.core)
        params_next = M_step_bi_fast(E_fit=E_next, params=params, X=x_k, Y=y_k, Z=Z, k=k)
        # params_next$B <- abs(params_next$B)
        
        if(i == 1) logl[1] = E_next$logl
        logl[i+1] = E_next$logl
        
        # # influence of Z
        # z_inf[i] <- with(params_next, sd(X%*%W/sig2E)/sd(Z%*%a/sig2G))
        
        # #debug
        # if(is.na(logl[i+1])) browser()
        # # if(logl[i+1] < logl[i]) browser()
        # if(logl[i+1] >1000) browser()
        fixtu[,k] <- E_next$mu_T[,k]
        fixtu[,r+k] <- E_next$mu_U[,k]
        
        if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
        if(i %in% c(1, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)) {
          print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
        }
        # solve(0)
        if(logl[i+1] > max(logl[1:i])) params_max <- params_next
        params = params_next
        
        
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
  
  outputt <- list(params = params) #, debug = debug_list)
  class(outputt) <- "PO2PLS"
  return(outputt)
}


Su_PO2PLS_bi_k <- function(X, Y, Z, r, rx, ry, k=1, steps = 1e5, tol = 1e-6, level=level, 
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
      E_next = E_step_bi_fast( X=x_k, Y=y_k, Z=Z, params=params, fixtu=fixtu, k=k, level=level, 
                               Nr.core=Nr.core)
      params_next = M_step_bi_fast(E_fit=E_next, params=params, X=x_k, Y=y_k, Z=Z, k=k)
      # params_next$B <- abs(params_next$B)
      
      if(i == 1) logl[1] = E_next$logl
      logl[i+1] = E_next$logl
      
      if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
      if(i %in% c(1, 1e1, 1e2, 1e3, 5e3, 1e4, 4e4)) {
        print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
      }
      # solve(0)
      if(logl[i+1] > max(logl[1:i])) params_max <- params_next
      params = params_next
    }
    if(!any(diff(logl[-1]) < -1e-10) | !random_restart_original) {
      random_restart = FALSE
      break
    }
    i_rr <- i_rr + 1
    params <- jitter_params(params)
    params[-(1:4)] <- generate_params_bi(X, Y, Z, r, rx, ry, type = 'r')[-(1:4)]
  }

  # latent variable values
  E_outp = E_step_bi_fast( X=x_k, Y=y_k, Z=Z, params=params, fixtu=fixtu, k=k, level=level, 
                  Nr.core=Nr.core)
  latent_var <- E_outp[names(E_outp) %in% c('mu_T', 'mu_U', 'mu_To', 'mu_Uo')]
  
  # R-square Z fit
  # mu_z <- with(params, a0 + X%*%W%*%(t(a)-B%*%t(b)) + Y%*%C%*%t(b))
  # z_p <- 1/(1+exp(-mu_z))
  # pre_z <- ifelse(z_p > 0.5, 1,0)
  # levs <- c(0,1)
  # zfit <- table(factor(pre_z, levs), factor(Z, levs))
  # acc <- caret::confusionMatrix(zfit)$overall['Accuracy']
  
  message("Nr steps was ", i)
  message("Log-likelihood: ", logl[i+1])
  # message("Accuracy of z fit: ", acc)
  outputt <- list(params = params, latent_var=latent_var, logl = logl[0:i+1][-1], 
                  GH_com=E_outp$GH_common,level=level) #, debug = debug_list)
  class(outputt) <- "PO2PLS"
  return(outputt)
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
  fixtu[,c(k,r+k)] <- mu_TU
  mu_TU <- fixtu
  
  # E(Sttuu|xyz) with numerical integration
  ## different part in integrand
  int_diff <- vector(mode = 'list', length = N)
  for(i in 1:N){
    int_diff[[i]] <- lapply(list_nodes, function(e){
      fixtu[i,c(k,r+k)] <- e
      crossprod(fixtu[i,,drop=F])
    })
  }
  ## Numerator
  S_ttuu_nu <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff[[e]], list_com_sub[[e]])))
  S_ttuu <- Map("/", S_ttuu_nu, dnmt)

  # Numerical estimation of Q_ab and gradient
  Q_ab <- GH_Q_ab_fast(cbind(a0,a,b), l_n=list_nodes, Z=Z, params=params, k=k, fixtu=fixtu, l_com=list_com_sub, de=dnmt)
  grd_ab <- GH_grd_ab_fast(l_n=list_nodes, Z=Z, params=params, k=k, fixtu=fixtu, l_com=list_com_sub, de=dnmt)
  # Backtracking rule to find the step size s
  s = 1
  Q_ab_new <- GH_Q_ab_fast((cbind(a0,a,b) + s*grd_ab),l_n=list_nodes, Z=Z, params=params, k=k, 
                           fixtu=fixtu, l_com=list_com_sub, de=dnmt)
  
  while(Q_ab_new < Q_ab + 0.5*s*tcrossprod(grd_ab)){
    s = 0.8*s
    # print(s)
    Q_ab_new <- GH_Q_ab_fast((cbind(a0,a,b) + s*grd_ab),l_n=list_nodes, Z=Z, params=params, k=k, 
                             fixtu=fixtu, l_com=list_com_sub, de=dnmt)
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

M_step_bi_fast <- function(E_fit, params, X, Y, Z, k){
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
    
    
    # # Stt, Suu different for each sample
    # tmp_sig2E = sapply(1:N, function(e){
    #   1/p * (sum(diag(X[e,,drop=F]%*%t(X[e,,drop=F]))) - 2*sum(diag(mu_T[e,,drop=F]%*%t(params$W)%*%t(X[e,,drop=F]))) +
    #            sum(diag(S_ttuu[[e]][1:r,1:r,drop=FALSE])) - sum(diag(params$SigTo)))
    # })
    # tmp_sig2F = sapply(1:N, function(e){
    #   1/q * (sum(diag(Y[e,,drop=F]%*%t(Y[e,,drop=F]))) - 2*sum(diag(mu_U[e,,drop=F]%*%t(params$C)%*%t(Y[e,,drop=F]))) +
    #            sum(diag(S_ttuu[[e]][r+1:r,r+1:r,drop=FALSE])) - sum(diag(params$SigUo)))
    # })
    
    # params$sig2E = mean(tmp_sig2E)
    # params$sig2F = mean(tmp_sig2F)
    
    # if(params$sig2E < 0) params$sig2E = 1e-5
    # if(params$sig2F < 0) params$sig2F = 1e-5
    
    params$B[k,k] = (t(Stu) %*% MASS::ginv(Stt))[k,k]
    params$SigT[k,k] = Stt[k,k]
    params$SigU[k,k] = Suu[k,k]
    params$SigH[k,k] = Shh[k,k]
    params$a0 = params$a0 + s*grd_ab[,1]
    params$a = params$a + s*grd_ab[,1+1:r]
    params$b = params$b + s*grd_ab[,1+r+1:r]
    
    params$W[,k] = orth(t(X/N) %*% mu_T[,k,drop=F] %*% (1/Stt[k,k,drop=F]))
    params$C[,k] = orth(t(Y/N) %*% mu_U[,k,drop=F] %*% (1/Suu[k,k,drop=F]))
    
    # params$sig2E = 1/p * abs(sum(diag(X%*%t(X)))/N - 2*sum(diag(mu_T%*%t(params$W)%*%t(X)))/N +
    #                        sum(diag(params$SigT)) - sum(diag(params$SigTo)))
    # params$sig2F = 1/q * abs(sum(diag(Y%*%t(Y)))/N - 2*sum(diag(mu_U%*%t(params$C)%*%t(Y)))/N +
    #                        sum(diag(params$SigU)) - sum(diag(params$SigUo)))
    
    return(params)
  })
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
  rr <- ncol(fixtu_i)/2
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


## function for numerical integration of Q_ab
GH_Q_ab_fast <- function(beta, l_n, Z, params, k, fixtu, l_com, de){
  N = length(Z)
  r=ncol(fixtu)/2
  # adjust (t,u) to (t,h)
  fixth <- fixtu
  fixth[,r+1:r] <- fixth[,r+1:r] - fixth[,1:r]%*% params$B
  l_n <- lapply(l_n, function(e) {
    e[,2] <- e[,2] - e[,1]%*% params$B[k,k,drop=F]
    return(e)
  })
  int_diff <- vector(mode = 'list', length = N)
  for(i in 1:N){
    int_diff[[i]] <- lapply(l_n, function(e){
      fixth[i,c(k,r+k)] <- e
      return(fixth[i,,drop=F])
    })
  }
  
  # log p(z|tu)
  int_diff_z <- vector(mode = 'list', length = N)
  for(i in 1:N){
    int_diff_z[[i]] <- lapply(int_diff[[i]], function(e){
      if(Z[i]==1){
        a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
        m <- max(0, a_tmp)
        -(m + log(exp(-m) + exp(a_tmp-m)))
      }else{
        a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
        m <- max(0, a_tmp)
        -(m + log(exp(-m) + exp(a_tmp-m)))
      }
    })
  }
  
  # combine common parts
  Q_ab <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff_z[[e]], l_com[[e]])))
  # Q_ab <- rapply(Q_ab, f=function(x) ifelse(is.nan(x),0,x), how="replace")
  Q_ab <- Map("/", Q_ab, de)
  # Add N samples
  Q_ab <- Reduce("+", Q_ab)
  return(Q_ab)
}

## function for numerical integration of gradient of Q_ab
GH_grd_ab_fast <- function(l_n, Z, params, k, fixtu, l_com, de){
  N = length(Z)
  r=ncol(fixtu)/2
  # adjust (t,u) to (t,h)
  fixth <- fixtu
  fixth[,r+1:r] <- fixth[,r+1:r] - fixth[,1:r]%*% params$B
  l_n <- lapply(l_n, function(e) {
    e[,2] <- e[,2] - e[,1]%*% params$B[k,k,drop=F]
    return(e)
  })
  int_diff <- vector(mode = 'list', length = N)
  for(i in 1:N){
    int_diff[[i]] <- lapply(l_n, function(e){
      fixth[i,c(k,r+k)] <- e
      return(fixth[i,,drop=F])
    })
  }
  
  beta <- with(params, cbind(a0,a,b))
  
  # log p(z|tu)
  int_diff_z <- vector(mode = 'list', length = N)
  for(i in 1:N){
    int_diff_z[[i]] <- lapply(int_diff[[i]], function(e){
      if(Z[i]==1){
        a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
        m <- max(0, a_tmp)
        exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)) * cbind(1,e)
      }else{
        a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
        m <- max(0, a_tmp)
        - exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)) * cbind(1,e)
      }
    })
  }
  
  # combine common parts
  grd_ab <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff_z[[e]], l_com[[e]])))
  # grd_ab <- rapply(grd_ab, f=function(x) ifelse(is.nan(x),0,x), how="replace")
  grd_ab <- Map("/", grd_ab, de)
  # Add N samples
  grd_ab <- Reduce("+", grd_ab)
  return(grd_ab)
}



# chi-squared test of (a0,a,b) for binary
chi2_ab_bi_fast <- function(fit, Z, k=1, level=NULL, Nr.core=NULL){
  B <- fit$params$B
  r = ncol(B)
  N = nrow(fit$latent_var$mu_T)
  fixtu <- with(fit$latent_var, cbind(mu_T, mu_U))
  
  # tmp.Estep <- E_step_bi(X, Y, Z, fit$par, level = level, Nr.core = Nr.core)
  # GH_com <- tmp.Estep$GH_common
  GH_com <- fit$GH_com
  
  list_com_log <- GH_com$list_com_log # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- GH_com$nodes
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
  
  # adjust (t,u) to (t,h)
  fixth <- fixtu
  fixth[,r+1:r] <- fixth[,r+1:r] - fixth[,1:r]%*% B
  list_nodes_th <- lapply(list_nodes, function(e) {
    e[,2] <- e[,2] - e[,1]%*% B[k,k,drop=F]
    return(e)
  })
  beta <- with(fit$params, cbind(a0,a,b))
  aabb <- with(fit$params, cbind(a,b))
  
  S_ab <- GH_I_ab_fast(beta,list_nodes_th,Z,list_com_sub,dnmt,'S',fixth,k)
  S2_ab <- GH_I_ab_fast(beta,list_nodes_th,Z,list_com_sub,dnmt,'S2',fixth,k)
  B_ab <- GH_I_ab_fast(beta,list_nodes_th,Z,list_com_sub,dnmt,'B',fixth,k)
  
  SiSj <- matrix(0, 1+2*r, 1+2*r)
  for(i in 2:N){
    for(j in 1:(i-1)){
      SiSj <- SiSj + S_ab[[i]]%*%t(S_ab[[j]])
    }
  }
  info_ab = (Reduce('+', B_ab) - Reduce('+', S2_ab) - SiSj - t(SiSj))[-1,-1]
  # print(info_ab)
  
  # test for all components
  stat = as.numeric(aabb %*% info_ab %*% t(aabb))
  if(stat<0){
    message('Negative Chi-square statistic for global test, try to set the level higher')
    p_global = NULL
  }else{
    p_global = pchisq(stat, df=2*r, lower.tail=FALSE)
  }
  
  # test for pair
  stat_pair <- p_pair <- c()
  for(k in 1:r){
    stat_pair[k] <- as.numeric(aabb[,c(k,k+r),drop=F] %*% info_ab[c(k,k+r),c(k,k+r),drop=F]
                            %*% t(aabb[,c(k,k+r),drop=F]))
    if(stat_pair[k]<0){
      stat_pair[k] = NULL
    } 
    else{
      p_pair[k] = pchisq(stat_pair[k], df=2, lower.tail=FALSE)
    }
  } 


  return(list(p_global=p_global, p_pair=p_pair))
}

## function for numerical integration of E[S(ab)|xyz], E[B(ab)|xyz], E[S(ab)S(ab)'|xyz]
GH_I_ab_fast <- function(beta, l_n=list_nodes_th,Z=Z,l_com=list_com_sub,de=dnmt,term,fixth,k){
  N = length(Z)
  r=ncol(fixth)/2

  int_diff <- vector(mode = 'list', length = N)
  for(i in 1:N){
    int_diff[[i]] <- lapply(l_n, function(e){
      fixth[i,c(k,r+k)] <- e
      return(fixth[i,,drop=F])
    })
  }
  
  int_diff_z <- vector(mode = 'list', length = N)
  
  if(term=="S"){  # Q'_ab
    for(i in 1:N){
      int_diff_z[[i]] <- lapply(int_diff[[i]], function(e){
        if(Z[i]==1){
          a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
          m <- max(0, a_tmp)
          exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)) * t(cbind(1,e))
        }else{
          a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
          m <- max(0, a_tmp)
          - exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)) * t(cbind(1,e))
        }
      })
    }
  }
  
  if(term=="B"){
    for(i in 1:N){
      int_diff_z[[i]] <- lapply(int_diff[[i]], function(e){
        if(Z[i]==1){
          a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
          m <- max(0, a_tmp)
          exp(a_tmp-m)/(exp(-m) + 2*exp(a_tmp-m) + exp(2*a_tmp-m)) * (t(cbind(1,e))%*%cbind(1,e))
        }else{
          a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
          m <- max(0, a_tmp)
          exp(a_tmp-m)/(exp(-m) + 2*exp(a_tmp-m) + exp(2*a_tmp-m)) * (t(cbind(1,e))%*%cbind(1,e))
        }
      })
    }
  }
  
  if(term=="S2"){
    for(i in 1:N){
      int_diff_z[[i]] <- lapply(int_diff[[i]], function(e){
        if(Z[i]==1){
          a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
          m <- max(0, a_tmp)
          (exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)))^2 * (t(cbind(1,e))%*%cbind(1,e))
        }else{
          a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
          m <- max(0, a_tmp)
          (exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)))^2 * (t(cbind(1,e))%*%cbind(1,e))
        }
      })
    }
  }
  

  # combine common parts
  int_ab <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff_z[[e]], l_com[[e]])))
  # grd_ab <- rapply(grd_ab, f=function(x) ifelse(is.nan(x),0,x), how="replace")
  int_ab <- Map("/", int_ab, de)
  return(int_ab)
}