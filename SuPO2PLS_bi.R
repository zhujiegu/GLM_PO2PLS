#' @export
E_step_bi <- function(X, Y, Z, params){
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
  a0 = a0
  a = params$a
  b = params$b
  ## define dimensions
  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  

  ## # E(tu|xyz) with numerical integration
  # N x r, collection of all the samples
  mu_TU <- lapply(1:N, function(e){
    GH_Intl(fun_mu, dim=2*r, level=6, X[e,],Y[e,],Z[e], params)}) %>% 
    unlist %>% matrix(nrow=2) %>% t
  
  # 2r x 2r, Average across N samples
  S_ttuu <- lapply(1:N, function(e){
    GH_Intl(fun_S, dim=2*r, level=6, X[e,],Y[e,],Z[e], params)})
  S_ttuu <- Reduce("+", S_ttuu)/N
  
  # loglikelihood
  lapply(1:N, function(e){
    GH_Intl(fun_com, dim=2*r, level=6, X[e,],Y[e,],Z[e], params)})
  
  list(
    mu_T = mu_T,
    mu_U = mu_U,
    mu_To = mu_To,
    mu_Uo = mu_Uo,
    Stt = Stt,
    Suu = Suu,
    Stoto = Stoto,
    Suouo = Suouo,
    Stu = Stu,
    Stto = Stto,
    Suuo = Suuo,
    Suto = Suto,
    See = Cee,
    Sff = Cff,
    Sgg = Cgg,
    Shh = Chh,
    aa=aa,
    bb=bb,
    EE=EE,
    ET=ET,
    loglik = loglik
  )
}

#' @export
M_step_su <- function(E_fit, params, X, Y, Z, orth_type = c("SVD","QR")){
  orth_x = ssq(params$Wo) > 0
  orth_y = ssq(params$Co) > 0
  orth_type = match.arg(orth_type)
  with(E_fit,{
    
    N = nrow(X)
    r = ncol(mu_T)
    rx = ncol(params$Wo)
    ry = ncol(params$Co)
    
    params$B = t(Stu) %*% MASS::ginv(Stt) * diag(1,r)
    params$SigT = Stt*diag(1,r)
    params$SigU = Suu*diag(1,r)
    params$SigTo = Stoto*diag(1,rx)
    params$SigUo = Suouo*diag(1,ry)
    params$SigH = Shh*diag(1,r)#abs(Suu - 2*Sut%*%params_old$B + Stt%*%params_old$B^2)
    
    params$sig2E = See
    params$sig2F = Sff
    params$sig2G = Sgg
    # params$sig2G = 0
    # params$sig2G = params$sig2G
    
    params$W = orth(t(X/N) %*% mu_T - params$Wo%*%t(Stto),type = orth_type)#%*%MASS::ginv(Stt)
    params$C = orth(t(Y/N) %*% mu_U - params$Co%*%t(Suuo),type = orth_type)#%*%MASS::ginv(Suu)
    
    params$Wo = suppressWarnings(orth_x*orth(t(X/N) %*% mu_To - params$W%*%Stto,type = orth_type))#%*%MASS::ginv(Stoto)
    params$Co = suppressWarnings(orth_y*orth(t(Y/N) %*% mu_Uo - params$C%*%Suuo,type = orth_type))#%*%MASS::ginv(Suouo)
    
    # params$a = (t(Z/N)%*%mu_T - params$b%*%Stu) %*% MASS::ginv(Stt*diag(1,r))
    # params$b = (t(Z/N)%*%mu_U - params$a%*%t(Stu)) %*% MASS::ginv(Suu*diag(1,r))
    
    params$a=aa
    params$b=bb
    
    # a + bB -> a
    # params$a = (t(Z/N)%*%mu_T - params$b%*%(Stu-Stt%*%B)) %*% MASS::ginv(Stt*diag(1,r))
    # params$b = (t(Z/N)%*%(mu_U-mu_T%*%B) - params$a%*%(Stu-Stt%*%B)) %*% MASS::ginv(Shh*diag(1,r))
    
    # A = (a,b)
    # Sttuu <- rbind(cbind(Stt*diag(1,r), Stu),
    #                cbind(t(Stu), Suu*diag(1,r)))
    # aabb <- t(Z/N) %*% cbind(mu_T,mu_U) %*% MASS::ginv(Sttuu)
    # params$a = aabb[,1:r, drop=FALSE]
    # params$b = aabb[,r+1:r, drop=FALSE]
    
    # # Try correct order
    # # browser()
    # signB <- sign(diag(params$B))
    # params$B <- params$B %*% diag(signB,r)
    # params$C <- params$C %*% diag(signB,r)
    # params$b <- params$b %*% diag(signB,r)
    # ordSB <- order(diag(params$SigT %*% params$B), decreasing = TRUE)
    # params$B <- params$B[,ordSB, drop=FALSE]
    # params$W <- params$W[,ordSB, drop=FALSE]
    # params$C <- params$C[,ordSB, drop=FALSE]
    # params$a <- params$a[,ordSB, drop=FALSE]
    # params$b <- params$b[,ordSB, drop=FALSE]
    # params$SigT <- params$SigT[,ordSB, drop=FALSE]
    # params$SigU <- params$SigU[,ordSB, drop=FALSE]
    # params$SigH <- params$SigH[,ordSB, drop=FALSE]
    
    # params$a = matrix(rep(0, r), nrow = 1, ncol=r)
    # params$b = matrix(rep(0, r), nrow = 1, ncol=r)
    
    return(params)
  })
}

#' @export
jitter_params <- function(params, amount = NULL){
  suppressWarnings(params[1:4] <- lapply(params[1:4], function(e) sign(ssq(e))*orth(jitter(e,amount = 1))))
  params
}


#' @export
Su_PO2PLS_bi <- function(X, Y, Z, r, rx, ry, steps = 1e5, tol = 1e-6, init_param= c('random','o2m'),
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
      params <- generate_params_su(X, Y, Z, r, rx, ry, type = "o2m")
    }
    else
    {
      if(r+max(rx,ry) > min(ncol(X),ncol(Y)) && init_param == "o2m")
      {
        cat("** NOTE: Too many components for init_param='o2m', switched to init_param='unit'**.\n\n");
        init_param = "unit"
      }
      params <- generate_params_su(X, Y, Z, r, rx, ry, type = init_param)
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
      E_next = E_step_bi(X, Y, Z, params)
      params_next = M_step_bi(E_next, params, X, Y, Z, orth_type = orth_type)
      params_next$B <- abs(params_next$B)
      
      if(i == 1) logl[1] = E_next$logl
      logl[i+1] = E_next$logl
      
      # # influence of Z
      # z_inf[i] <- with(params_next, sd(X%*%W/sig2E)/sd(Z%*%a/sig2G))
      
      # #debug
      # if(is.na(logl[i+1])) browser()
      # # if(logl[i+1] < logl[i]) browser()
      
      if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
      if(i %in% c(1e1, 1e2, 1e3, 5e3, 1e4, 4e4)) {
        print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
      }
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
  signB <- sign(diag(params$B))
  params$B <- params$B %*% diag(signB,r)
  params$C <- params$C %*% diag(signB,r)
  params$b <- params$b %*% diag(signB,r)
  
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