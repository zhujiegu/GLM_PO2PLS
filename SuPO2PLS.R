
#' @export
blockm<-function(A,B,C)
  #input: Matrices A,B,C
  #output: the block matrix
  # A    B
  #t(B)  C
{
  M = rbind(cbind(A,B),cbind(t(B),C))
  return(M)
}

#' @export
generate_params <- function(X, Y, Z, r, rx, ry, alpha = 0.1, type=c('o2m','random')){
  type=match.arg(type)
  p = ifelse(is.matrix(X) | type != "random", ncol(X), X)
  q = ifelse(is.matrix(Y) | type != "random", ncol(Y), Y)
  if(type=="o2m"){
    fit <- o2m(X, Y, r, rx, ry, stripped=TRUE)
    return(with(fit, {
      x_tp <- cbind(Tt,U)
      ab <- MASS::ginv((t(x_tp) %*% x_tp)) %*% t(x_tp) %*% Z
      a = t(ab[1:r,])
      b = t(ab[(r+1):(2*r),])
      list(
        W = W.,
        Wo = suppressWarnings(orth(P_Yosc.)),
        C = C.,
        Co = suppressWarnings(orth(P_Xosc.)),
        B = abs(cov(Tt,U)%*%MASS::ginv(cov(Tt)))*diag(1,r),
        SigT = cov(Tt)*diag(1,r),
        SigTo = sign(rx)*cov(T_Yosc.)*diag(1,max(1,rx)),
        SigUo = sign(ry)*cov(U_Xosc.)*diag(1,max(1,ry)),
        SigH = cov(H_UT)*diag(1,r),
        sig2E = (ssq(X)-ssq(Tt)-ssq(T_Yosc.))/prod(dim(X)) + 0.01,
        sig2F = (ssq(Y)-ssq(U)-ssq(U_Xosc.))/prod(dim(Y)) + 0.01,
        a = a,
        b = b,
        sig2G = as.numeric(cov(Z-Tt%*%t(a)-U%*%t(b)))
      )}))
  }
  if(type=="random"){
    outp <- list(
      W = orth(matrix(rnorm(p*r), p, r)+1),
      Wo = suppressWarnings(sign(rx)*orth(matrix(rnorm(p*max(1,rx)), p, max(1,rx))+seq(-p/2,p/2,length.out = p))),
      C = orth(matrix(rnorm(q*r), q, r)+1),
      Co = suppressWarnings(sign(ry)*orth(matrix(rnorm(q*max(1,rx)), q, max(1,ry))+seq(-q/2,q/2,length.out = q))),
      B = diag(sort(runif(r,1,4),decreasing = TRUE),r),
      SigT = diag(sort(runif(r,1,3),decreasing = TRUE),r),
      SigTo = sign(rx)*diag(sort(runif(max(1,rx),1,3),decreasing = TRUE),max(1,rx)),
      SigUo = sign(ry)*diag(sort(runif(max(1,ry),1,3),decreasing = TRUE),max(1,ry)),
      a = matrix(rnorm(r), nrow = 1, ncol = r),
      b = matrix(rnorm(r), nrow = 1, ncol = r)
    )
    outp$SigH = diag(alpha/(1-alpha)*(diag(outp$SigT%*%outp$B^2)), r) #cov(H_UT)*diag(1,r),
    outp$SigU <- with(outp, SigT %*% B^2 + SigH)
    return(with(outp, {
      c(outp,
        sig2E = alpha/(1-alpha)*(mean(diag(SigT)) + mean(diag(SigTo)))/p,
        sig2F = alpha/(1-alpha)*(mean(diag(SigT%*%B^2 + SigH)) + mean(diag(SigUo)))/q,
        sig2G = as.numeric(alpha/(1-alpha)*(a %*% SigT %*% t(a) + b %*% SigU %*% t(b))))
    }))
  }
}

#' @export
generate_data <- function(N, params, distr = rnorm){
  W = params$W
  C = params$C
  Wo = params$Wo
  Co = params$Co
  B = params$B
  SigT = params$SigT
  SigTo = params$SigTo + 1e-6*SigT[1]*(params$SigTo[1]==0)
  SigH = params$SigH
  sig2E = params$sig2E
  sig2F = params$sig2F
  SigU = SigT%*%B^2 + SigH
  SigUo = params$SigUo + 1e-6*SigU[1]*(params$SigUo[1]==0)
  a = params$a
  b = params$b
  sig2G = params$sig2G
  
  
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co),
                cbind(a, b, matrix(0,1,rx+ry)))
  VarM = blockm(
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)
  
  # MASS::mvrnorm(n = N,
  #               mu = rep(0,p+q),
  #               Sigma = Gamma %*% VarZ %*% t(Gamma) +
  #                 diag(rep(c(sig2E,sig2F),c(p,q))))
  
  M <- scale(matrix(rnorm(N*(2*r+rx+ry)), N))
  M <- M %*% chol(VarM)
  M[,2*r+1:rx] <- sign(ssq(Wo))*M[,2*r+1:rx]
  M[,2*r+rx+1:ry] <- sign(ssq(Co))*M[,2*r+rx+1:ry]
  
  EFG <- cbind(scale(matrix(rnorm(N*p), N))*sqrt(sig2E), 
               scale(matrix(rnorm(N*q), N))*sqrt(sig2F),
               scale(matrix(rnorm(N), N))*sqrt(sig2G))
  
  M %*% t(Gamma) + EFG
  
}

#' @export
Lemma <- function(XYZ, Mtilde, GammaEFG, p, q, r, rx, ry){
  
  EE <- XYZ %*% (GammaEFG %*% Mtilde)
  ET <- Mtilde
  mu_T <- EE[,1:r]
  mu_U <- EE[,(r+1):(2*r)]
  mu_To <- EE[,(2*r+1):(2*r+rx)]
  mu_Uo <- EE[,(2*r+rx+1):(2*r+rx+ry)]
  Stt <- ET[1:r,1:r]
  Suu <- ET[(r+1):(2*r), (r+1):(2*r)]
  Stoto <- ET[(2*r+1):(2*r+rx), (2*r+1):(2*r+rx)]
  Suouo <- ET[(2*r+rx+1):(2*r+rx+ry), (2*r+rx+1):(2*r+rx+ry)]
  Stu <- ET[1:r, (r+1):(2*r)]
  Stto <- ET[1:r, (2*r+1):(2*r+rx)]
  Suto <- ET[(r+1):(2*r), (2*r+1):(2*r+rx)]
  Suuo <- ET[(r+1):(2*r), (2*r+rx+1):(2*r+rx+ry)]
  return(list(mu_T=mu_T,mu_U=mu_U,mu_To=mu_To,mu_Uo=mu_Uo,Stt=Stt,Suu=Suu,Stoto=Stoto,
              Suouo=Suouo,Stu=Stu,Stto=Stto,Suto=Suto,Suuo=Suuo))
}


#' @export
E_step <- function(X, Y, Z, params){
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
  a = params$a
  b = params$b
  sig2G = params$sig2G

  ## define dimensions
  N = nrow(X)
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  
  ## concatenate data
  dataXYZ <- cbind(X,Y,Z)
  
  ## Gamma is the generalized loading matrix, with PO2PLS structure
  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co),
                cbind(a, b, matrix(0,1,rx+ry)))
  ## Gamma multiplied by inverse SigmaEF
  GammaEFG <- Gamma
  GammaEFG[1:p,c(1:r,2*r+1:rx)] <- MASS::ginv(sig2E)[1]* GammaEFG[1:p,c(1:r,2*r+1:rx)]
  GammaEFG[p+1:q,c(r+1:r,2*r+rx+1:ry)] <- MASS::ginv(sig2F)[1]* GammaEFG[p+1:q,c(r+1:r,2*r+rx+1:ry)]
  GammaEFG[p+q+1,1:(2*r)] <- MASS::ginv(sig2G)[1]* GammaEFG[p+q+1,1:(2*r)]
  # all.equal(GammaEFG, diag(c(rep(1/sig2E,p),rep(1/sig2F,q),1/sig2G))%*%Gamma)
  GGef <- t(Gamma) %*% GammaEFG
  
  
  ## diagonal cov matrix of (E,F), hopefully NOT NEEDED
  # SigmaEF = diag(rep(c(sig2E,sig2F),c(p,q)))
  ## ALMOST diagonal cov matrix of (T,U,To,Uo)
  SigmaM = blockm(
    blockm(
      blockm(SigT, SigT%*%B, SigU),
      matrix(0,2*r,rx), SigTo),
    matrix(0,2*r+rx,ry), SigUo)
  
  ## inverse middle term lemma
  Mtilde <- MASS::ginv(MASS::ginv(SigmaM) + GGef)

  ## Calculate conditional expectations with efficient lemma
  EE <- dataXYZ %*% (GammaEFG %*% Mtilde)
  ET <- Mtilde
  # print(Mtilde)
  # print(GGef)
  # print(SigmaM)
  mu_T <- EE[,1:r,drop=FALSE]
  mu_U <- EE[,(r+1):(2*r),drop=FALSE]
  mu_To <- EE[,(2*r+1):(2*r+rx),drop=FALSE]
  mu_Uo <- EE[,(2*r+rx+1):(2*r+rx+ry),drop=FALSE]
  Stt <- ET[1:r,1:r,drop=FALSE] + crossprod(mu_T)/N
  Suu <- ET[(r+1):(2*r), (r+1):(2*r),drop=FALSE] + crossprod(mu_U)/N
  Stoto <- ET[(2*r+1):(2*r+rx), (2*r+1):(2*r+rx),drop=FALSE] + crossprod(mu_To)/N
  Suouo <- ET[(2*r+rx+1):(2*r+rx+ry), (2*r+rx+1):(2*r+rx+ry),drop=FALSE] + crossprod(mu_Uo)/N
  Stu <- ET[1:r, (r+1):(2*r),drop=FALSE] + crossprod(mu_T,mu_U)/N
  Stto <- ET[1:r, (2*r+1):(2*r+rx),drop=FALSE] + crossprod(mu_T, mu_To)/N
  Suto <- ET[(r+1):(2*r), (2*r+1):(2*r+rx),drop=FALSE] + crossprod(mu_U, mu_To)/N
  Suuo <- ET[(r+1):(2*r), (2*r+rx+1):(2*r+rx+ry),drop=FALSE] + crossprod(mu_U, mu_Uo)/N
  
  ## Calculate cond mean of E,F,G,H
  # invS_covEF <- diag(1,p+q) - invEF_Gamma %*% invZtilde %*% t(Gamma)
  # covEF = rbind(diag(sig2E,p), diag(0,q,p))
  # mu_EF_old = dataXY - dataXY %*% invEF_Gamma %*% invZtilde %*% t(Gamma)
  mu_EFG = dataXYZ
  mu_EFG <- mu_EFG - (dataXYZ %*% (GammaEFG %*% Mtilde)) %*% t(Gamma)
  
  ## Calculate immediately expected crossprod of E,F
  
  # Cee_old <- sum(diag(Ceeff_old[1:p,1:p]))/p
  Cee <- sum(diag(
    crossprod(rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                    cbind(matrix(0,q,r), 0*C, matrix(0,q,rx), 0*Co)))%*%Mtilde
  ))/p + ssq(mu_EFG[,1:p])/N/p
  
  # Slow, directly from the formula
  # Cee_test <- t(X)%*%X/N - 2*t(X)%*%mu_T%*%t(W)/N -2*t(X)%*%mu_To%*%t(Wo)/N + W%*%Stt%*%t(W) + 2*W%*%Stto%*%t(Wo) + Wo%*%Stoto%*%t(Wo)
  # Cee_test <- (Cee_test %>% diag %>% sum)/p
  # all.equal(Cee, Cee_test)
  
  Cff <- sum(diag(
    crossprod(rbind(cbind(0*W, matrix(0,p,r), 0*Wo, matrix(0,p,ry)),
                    cbind(matrix(0,q,r), C, matrix(0,q,rx), Co)))%*%Mtilde
  ))/q + ssq(mu_EFG[,p+1:q])/N/q
  
  # Slow, directly from the formula
  # Cff_test <- t(Y)%*%Y/N - 2*t(Y)%*%mu_U%*%t(C)/N -2*t(Y)%*%mu_Uo%*%t(Co)/N + C%*%Suu%*%t(C) + 2*C%*%Suuo%*%t(Co) + Co%*%Suouo%*%t(Co)
  # Cff_test <- (Cff_test %>% diag %>% sum)/q
  # all.equal(Cff, Cff_test)
  
  Chh <- Suu - t(Stu)%*%B - t(B)%*%Stu + t(B)%*%Stt%*%B
  
  # covH = rbind(0*W, C%*%SigH, b%*%SigH)
  # covHEFG = rbind(0*W, C%*%SigH/sig2F, b%*%SigH/sig2G)
  # # invS_covH <- (covH/sig2F - invEF_Gamma %*% invZtilde %*% crossprod(invEF_Gamma,covH))
  # # mu_H_old = dataXY %*% invS_covH
  # # mu_H_old <- dataXY %*% invS %*% covH
  # mu_H <- dataXYZ %*% covHEFG
  # mu_H <- mu_H - (dataXYZ %*% (GammaEFG %*% Mtilde)) %*% (t(Gamma) %*% covHEFG)
  # # Chh_old = SigH - t(covH) %*% invS_covH + crossprod(mu_H) / N
  # # Chh_old <- SigH - t(covH) %*% invS %*% covH + crossprod(mu_H) / N
  # Chh <- SigH
  # Chh <- Chh - t(covH) %*% covHEFG
  # Chh <- Chh + (t(covH) %*% GammaEFG %*% Mtilde) %*% (t(Gamma) %*% covHEFG)
  # Chh <- Chh + crossprod(mu_H) / N

  # (a,b) -> A
  Sttuu <- rbind(cbind(Stt, Stu),
                 cbind(t(Stu), Suu))
  aabb <- t(Z/N) %*% cbind(mu_T,mu_U) %*% MASS::ginv(Sttuu)
  
  aa = aabb[,1:r, drop=FALSE]
  bb = aabb[,r+1:r, drop=FALSE]
  
  Cgg <- as.numeric(crossprod(Z)/N - 2*t(Z)%*%cbind(mu_T,mu_U)%*%t(aabb)/N +
                          aabb%*%Sttuu%*%t(aabb))
  
  # print(all.equal(Cgg,Cgg_old))
  # a+bB -> a
  # Cgg <- as.numeric(crossprod(Z)/N - 2*t(Z)%*%mu_T%*%t(a)/N -2*t(Z)%*%(mu_U-mu_T%*%B)%*%t(b)/N +
  #                     a%*%Stt%*%t(a) + 2*a%*%(Stu-Stt%*%B)%*%t(b) +
  #                     b%*%(Suu - t(Stu)%*%B - t(B)%*%Stu + t(B)%*%Stt%*%B)%*%t(b))
  
  ## log of det SigmaXY, see matrix determinant lemma
  logdet <- log(det(diag(2*r+rx+ry) + GGef%*%SigmaM))+p*log(sig2E)+q*log(sig2F)
  logdet <- logdet + log(ifelse(sig2G>0, sig2G, 1)) 
  XYZinvS <- ssq(cbind(X/sqrt(sig2E), Y/sqrt(sig2F), Z*MASS::ginv(sqrt(sig2G))[1]))
  XYZinvS <- XYZinvS - sum(diag(crossprod(dataXYZ %*% GammaEFG) %*% Mtilde))
  ## Log likelihood
  loglik = N*(p+q+1)*log(2*pi) + N * logdet + XYZinvS
  loglik = - loglik/2
  
  # # Direct log-likelihood
  # Sig <- rbind(cbind(W%*%SigT%*%t(W) + Wo%*%SigTo%*%t(Wo) + sig2E*diag(p), W%*%SigT%*%B%*%t(C), W%*%SigT%*%(t(a) + B%*%t(b))),
  #              cbind(C%*%B%*%SigT%*%t(W), C%*%SigU%*%t(C) + Co%*%SigUo%*%t(Co) + sig2F*diag(q), C%*%SigT%*%B%*%t(a) + C%*%SigU%*%t(b)),
  #              cbind((a+b%*%B)%*%SigT%*%t(W), a%*%B%*%SigT%*%t(C) + b%*%SigU%*%t(C), a%*%SigT%*%t(a) + b%*%SigU%*%t(b) + 2*a%*%B%*%SigT%*%t(b) + sig2G))
  # 
  # invSig <- MASS::ginv(Sig)
  # loglik_test = N*(p+q+1)*log(2*pi) + N*log(det(Sig)) + sum(diag((crossprod(dataXYZ)%*%invSig)))
  # loglik_test = -loglik_test/2
  # print(all.equal(loglik, loglik_test))
  
  # Direct log-likelihood a+bB -> a
  # Sig <- rbind(cbind(W%*%SigT%*%t(W) + Wo%*%SigTo%*%t(Wo) + sig2E*diag(p), W%*%SigT%*%B%*%t(C), W%*%SigT%*%t(a)),
  #              cbind(C%*%B%*%SigT%*%t(W), C%*%SigU%*%t(C) + Co%*%SigUo%*%t(Co) + sig2F*diag(q), C%*%(SigH%*%t(b) + B%*%SigT%*%t(a))),
  #              cbind(a%*%SigT%*%t(W), (a%*%SigT%*%B+b%*%SigH)%*%t(C), a%*%SigT%*%t(a) + b%*%SigH%*%t(b) + sig2G))
  # 
  # invSig <- MASS::ginv(Sig)
  # loglik_test = N*(p+q+1)*log(2*pi) + N*log(det(Sig)) + sum(diag((crossprod(dataXYZ)%*%invSig)))
  # loglik_test = -loglik_test/2
  # 
  # print(loglik)
  # print(loglik_test)
  # all.equal(loglik, loglik_test)

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
    loglik = loglik,
    EE = EE,
    GammaEFG = GammaEFG,
    Mtilde=Mtilde
  )
}

#' @export
M_step <- function(E_fit, params, X, Y, Z, orth_type = c("SVD","QR")){
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
diagnostics.PO2PLS <- function(th, th0){
  c(
    W = max(abs(crossprod(th$W,th0$W))),
    C = max(abs(crossprod(th$C,th0$C))),
    Wo = max(abs(crossprod(th$Wo,th0$Wo))),
    Co = max(abs(crossprod(th$Co,th0$Co))),
    varTo_T = sum(diag(th$SigTo))/sum(diag(th$SigT))/ncol(th$Wo)*ncol(th$W),
    varUo_U = sum(diag(th$SigUo))/sum(diag(th$SigT%*%th$B+th$SigH))/ncol(th$Co)*ncol(th$C),
    varU_T = sum(diag(th$SigT%*%th$B+th$SigH))/sum(diag(th$SigT))
  )
}

#' @export
Su_PO2PLS <- function(X, Y, Z, r, rx, ry, steps = 1e5, tol = 1e-6, init_param= c('random','o2m'),
                   orth_type = "SVD", random_restart = FALSE){
  # if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F") %in% names(init_param))) {cat('using old fit \n'); params <- init_param}
  # else {params <- generate_params(X, Y, Z, r, rx, ry, type = init_param)}
  
  debug_list <- list()
  if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a","b","sig2G") %in% names(init_param))){
    message('using old fit \n')
    params <- init_param}
  else {
    init_param <- match.arg(init_param)
    if(r+max(rx,ry) <= min(ncol(X),ncol(Y)) && init_param == "o2m")
    {
      params <- generate_params(X, Y, Z, r, rx, ry, type = "o2m")
    }
    else
    {
      if(r+max(rx,ry) > min(ncol(X),ncol(Y)) && init_param == "o2m")
      {
        cat("** NOTE: Too many components for init_param='o2m', switched to init_param='unit'**.\n\n");
        init_param = "unit"
      }
      params <- generate_params(X, Y, Z, r, rx, ry, type = init_param)
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
      E_next = E_step(X, Y, Z, params)
      params_next = M_step(E_next, params, X, Y, Z, orth_type = orth_type)
      params_next$B <- abs(params_next$B)

      if(i == 1) logl[1] = E_next$logl
      logl[i+1] = E_next$logl
      
      #debug
      if(is.na(logl[i+1])) browser()
      # if(logl[i+1] < logl[i]) browser()
      
      if(i > 1 && abs(logl[i+1]-logl[i]) < tol) break
      if(i %in% c(1e1, 1e2, 1e3, 5e3, 1e4, 4e4)) {
        print(data.frame(row.names = 1, steps = i, time = unname(proc.time()-tic)[3], diff = logl[i+1]-logl[i], logl = logl[i+1]))
      }
      if(logl[i+1] > max(logl[1:i])) params_max <- params_next
      params = params_next
      
      debug_list[[i]] <- list(Cor=cor(cbind(E_next$mu_T,E_next$mu_U)),
                              aa=params_next$a,
                              bb=params_next$b,
                              sig2G=params_next$sig2G,
                              Sttuu=with(E_next, rbind(cbind(2*Stt*diag(1,r), Stu+t(Stu)),
                                             cbind(Stu+t(Stu), 2*Suu*diag(1,r)))),
                              ZTU=t(Z/N)%*%cbind(E_next$mu_T,E_next$mu_U))
    }
    if(!any(diff(logl[-1]) < -1e-10) | !random_restart_original) {
      random_restart = FALSE
      break
    }
    i_rr <- i_rr + 1
    params <- jitter_params(params)
    params[-(1:4)] <- generate_params(X, Y, r, rx, ry, type = 'r')[-(1:4)]
  }
  # params <- params_max
  signB <- sign(diag(params$B))
  params$B <- params$B %*% diag(signB,r)
  params$C <- params$C %*% diag(signB,r)
  params$b <- params$b %*% diag(signB,r)
  ordSB <- order(diag(params$SigT %*% params$B), decreasing = TRUE)
  params$B <- params$B[ordSB,ordSB, drop=FALSE]
  params$W <- params$W[,ordSB, drop=FALSE]
  params$C <- params$C[,ordSB, drop=FALSE]
  params$a <- params$a[,ordSB, drop=FALSE]
  params$b <- params$b[,ordSB, drop=FALSE]
  params$SigT <- params$SigT[ordSB,ordSB, drop=FALSE]
  params$SigU <- params$SigU[ordSB,ordSB, drop=FALSE]
  params$SigH <- params$SigH[ordSB,ordSB, drop=FALSE]
  
  message("Nr steps was ", i)
  message("Negative increments: ", any(diff(logl[1:i+1]) < 0), "\n",
          "Smallest increments: ", min(diff(logl[1:i+1])),
          "; Last increment: ", signif(logl[i+1]-logl[i],4))
  message("Log-likelihood: ", logl[i+1])
  outputt <- list(params = params, logl = logl[0:i+1][-1], debug = debug_list)
  class(outputt) <- "PO2PLS"
  return(outputt)
}
