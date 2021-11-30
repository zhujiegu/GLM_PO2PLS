
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
#' Provide XYZ if "o2m". Provide p,q if "random" or "specify"
#' Note in "specify", set.seed is used
generate_params_su <- function(X=NULL, Y=NULL, Z=NULL, p=NULL, q=NULL, r, rx, ry, alpha_x = 0.1, alpha_y = 0.1, alpha_z = 0.1,
                               alpha_tu = 0.1, B=1, a=t(2),b=t(1), type=c('o2m','random','specify')){
  type=match.arg(type)
  if(type=="o2m"){
    p = ncol(X)
    q = ncol(Y)
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
        SigTo = sign(rx)*cov(T_Yosc)*diag(1,max(1,rx)),
        SigUo = sign(ry)*cov(U_Xosc)*diag(1,max(1,ry)),
        SigH = cov(H_UT)*diag(1,r),
        sig2E = (ssq(X)-ssq(Tt)-ssq(T_Yosc))/prod(dim(X)) + 0.01,
        sig2F = (ssq(Y)-ssq(U)-ssq(U_Xosc))/prod(dim(Y)) + 0.01,
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
      B = diag(B, r),
      # B = diag(sort(runif(r,1,1),decreasing = TRUE),r), # set to 1 for simulation
      SigT = diag(sort(runif(r,1,3),decreasing = TRUE),r),
      SigTo = sign(rx)*diag(sort(runif(max(1,rx),1,3),decreasing = TRUE),max(1,rx)),
      SigUo = sign(ry)*diag(sort(runif(max(1,ry),1,3),decreasing = TRUE),max(1,ry)),
      a = a,
      b = b
    )
    outp$SigH = diag(alpha_tu/(1-alpha_tu)*(diag(outp$SigT%*%outp$B^2)), r) #cov(H_UT)*diag(1,r),
    outp$SigU <- with(outp, SigT %*% B^2 + SigH)
    return(with(outp, {
      c(outp,
        sig2E = alpha_x/(1-alpha_x)*(mean(diag(SigT)) + mean(diag(SigTo)))/p,
        sig2F = alpha_y/(1-alpha_y)*(mean(diag(SigT%*%B^2 + SigH)) + mean(diag(SigUo)))/q,
        sig2G = as.numeric(alpha_z/(1-alpha_z)*(a %*% SigT %*% t(a) + b %*% SigU %*% t(b))))
    }))
  }
  if(type=="specify"){
    set.seed(0)
    outp <- list(
      W = orth(matrix(rnorm(p*r), p, r)+1),
      Wo = suppressWarnings(sign(rx)*orth(matrix(rnorm(p*max(1,rx)), p, max(1,rx))+seq(-p/2,p/2,length.out = p))),
      C = orth(matrix(rnorm(q*r), q, r)+1),
      Co = suppressWarnings(sign(ry)*orth(matrix(rnorm(q*max(1,rx)), q, max(1,ry))+seq(-q/2,q/2,length.out = q))),
      B = diag(B, r),
      # B = diag(sort(runif(r,1,1),decreasing = TRUE),r), # set to 1 for simulation
      SigT = diag(3,r),
      SigTo = sign(rx)*diag(3,max(1,rx)),
      SigUo = sign(ry)*diag(3,max(1,ry)),
      a = a,
      b = b
    )
    outp$SigH = diag(alpha_tu/(1-alpha_tu)*(diag(outp$SigT%*%outp$B^2)), r) #cov(H_UT)*diag(1,r),
    outp$SigU <- with(outp, SigT %*% B^2 + SigH)
    return(with(outp, {
      c(outp,
        sig2E = alpha_x/(1-alpha_x)*(mean(diag(SigT)) + mean(diag(SigTo)))/p,
        sig2F = alpha_y/(1-alpha_y)*(mean(diag(SigT%*%B^2 + SigH)) + mean(diag(SigUo)))/q,
        sig2G = as.numeric(alpha_z/(1-alpha_z)*(a %*% SigT %*% t(a) + b %*% SigU %*% t(b))))
    }))
  }
}

#' @export
generate_data_su <- function(N, params, distr = rnorm){
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
                cbind(a-b%*%B, b, matrix(0,1,rx+ry)))
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
Lemma_su <- function(XYZ, Mtilde, GammaEFG, p, q, r, rx, ry){
  
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
E_step_su <- function(X, Y, Z, params, b_on_h=F){
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
  if(!b_on_h){
    Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                  cbind(matrix(0,q,r), C, matrix(0,q,rx), Co),
                  cbind(a, b, matrix(0,1,rx+ry)))
  }else{
    Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                  cbind(matrix(0,q,r), C, matrix(0,q,rx), Co),
                  cbind(a-b%*%B, b, matrix(0,1,rx+ry)))
  }

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

  # (a,b) -> A (correlated TU)
  if(!b_on_h){
    Sttuu <- rbind(cbind(Stt, Stu),
                   cbind(t(Stu), Suu))
    aabb <- t(Z/N) %*% cbind(mu_T,mu_U) %*% MASS::ginv(Sttuu)
    
    aa = aabb[,1:r, drop=FALSE]
    bb = aabb[,r+1:r, drop=FALSE]
    
    Cgg <- as.numeric(crossprod(Z)/N - 2*t(Z)%*%cbind(mu_T,mu_U)%*%t(aabb)/N +
                        aabb%*%Sttuu%*%t(aabb)) # for some reason, we need to use updated ab here
  }else{   # Uncorrelate TU
    Sth <- Stu - Stt%*%B
    mu_H = mu_U - mu_T %*% B
    Stthh <- rbind(cbind(Stt, Sth),
                   cbind(t(Sth), Chh))
    aabb <- t(Z/N) %*% cbind(mu_T,mu_H) %*% MASS::ginv(Stthh)
    
    aa = aabb[,1:r, drop=FALSE]
    bb = aabb[,r+1:r, drop=FALSE]
    
    # Cgg <- sum(diag(
    #   crossprod(cbind(a-b%*%B, b, matrix(0,1,rx+ry)))%*%Mtilde
    # )) + ssq(mu_EFG[,p+q+1])/N
    Cgg <- as.numeric(crossprod(Z)/N - 2*t(Z)%*%cbind(mu_T,mu_H)%*%t(cbind(a,b))/N +
                        cbind(a,b)%*%Stthh%*%t(cbind(a,b)))
  }

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
  
  # # Direct log-likelihood b_on_h
  # Sig <- rbind(cbind(W%*%SigT%*%t(W) + Wo%*%SigTo%*%t(Wo) + sig2E*diag(p), W%*%SigT%*%B%*%t(C), W%*%SigT%*%(t(a))),
  #              cbind(C%*%B%*%SigT%*%t(W), C%*%SigU%*%t(C) + Co%*%SigUo%*%t(Co) + sig2F*diag(q), C%*%B%*%SigT%*%t(a) + C%*%SigH%*%t(b)),
  #              cbind(a%*%SigT%*%t(W), a%*%SigT%*%B%*%t(C) + b%*%SigH%*%t(C), a%*%SigT%*%t(a) + b%*%SigH%*%t(b) + sig2G))
  # invSig <- MASS::ginv(Sig)
  # loglik_test = N*(p+q+1)*log(2*pi) + N*log(det(Sig)) + sum(diag((crossprod(dataXYZ)%*%invSig)))
  # loglik_test = -loglik_test/2
  # print(loglik)
  # print(loglik_test)
  # print(all.equal(loglik, loglik_test))
  
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
Su_PO2PLS <- function(X, Y, Z, r, rx, ry, steps = 1e5, tol = 1e-6, init_param= c('random','o2m'),
                   orth_type = "SVD", random_restart = FALSE, b_on_h=T){
  # if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F") %in% names(init_param))) {cat('using old fit \n'); params <- init_param}
  # else {params <- generate_params(X, Y, Z, r, rx, ry, type = init_param)}
  
  # debug_list <- list()
  # z_inf <- c()
  if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a","b","sig2G") %in% names(init_param))){
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
      E_next = E_step_su(X, Y, Z, params, b_on_h = b_on_h)
      params_next = M_step_su(E_next, params, X, Y, Z, orth_type = orth_type)
      params_next$B <- abs(params_next$B)

      if(i == 1) logl[1] = E_next$logl
      logl[i+1] = E_next$logl
      
      # # influence of Z
      # z_inf[i] <- with(params_next, sd(X%*%W/sig2E)/sd(Z%*%a/sig2G))
      
      #debug
      if(is.na(logl[i+1])) browser()
      # if(logl[i+1] < logl[i]) browser()
      
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
  ordSB <- order(diag(params$SigT %*% params$B), decreasing = TRUE)
  params$B <- params$B[ordSB,ordSB, drop=FALSE]
  params$W <- params$W[,ordSB, drop=FALSE]
  params$C <- params$C[,ordSB, drop=FALSE]
  params$a <- params$a[,ordSB, drop=FALSE]
  params$b <- params$b[,ordSB, drop=FALSE]
  params$SigT <- params$SigT[ordSB,ordSB, drop=FALSE]
  params$SigU <- params$SigU[ordSB,ordSB, drop=FALSE]
  params$SigH <- params$SigH[ordSB,ordSB, drop=FALSE]
  
  # R-square Z fit 
  pre_z <- with(params, X%*%W%*%(t(a)-B%*%t(b)) + Y%*%C%*%t(b))
  zfit <- as.numeric(cor(pre_z,Z)^2)
  
  message("Nr steps was ", i)
  message("Negative increments: ", any(diff(logl[1:i+1]) < 0), "\n",
          "Smallest increments: ", min(diff(logl[1:i+1])),
          "; Last increment: ", signif(logl[i+1]-logl[i],4))
  message("Log-likelihood: ", logl[i+1])
  message("R-square of z fit: ", zfit)
  outputt <- list(params = params, logl = logl[0:i+1][-1], zfit = zfit) #, debug = debug_list)
  class(outputt) <- "PO2PLS"
  return(outputt)
}

sd_B <- function(fit, X, Y){
  # tmp.Estep <- E_step(X, Y, fit$par)
  # with(tmp.Estep,
  #      Stt%*%solve(fit$par$SigH) -
  #        (crossprod(Sut) - crossprod(Stt)%*%fit$par$B^2)%*%
  #        solve(fit$par$SigH^2)) %>%
  #   multiply_by(nrow(X)) %>% solve %>% diag %>% abs %>% raise_to_power(0.5)
  
  tmp.Estep <- E_step_su(X, Y, Z, fit$par)
  with(tmp.Estep,
       solve(fit$par$SigH)%*%Stt - solve(fit$par$SigH)%*%
         (Sut%*%t(Sut) - Sut%*%Stt%*%fit$par$B - t(Sut%*%Stt%*%fit$par$B)
          + fit$par$B%*%crossprod(Stt)%*%fit$par$B) %*%solve(fit$par$SigH)) %>%
    multiply_by(nrow(X)) %>% solve %>% diag %>% abs %>% raise_to_power(0.5)
}

sd_ab <- function(fit, X, Y, Z, b_on_h = T){
  if(!b_on_h) stop("please set b_on_h to TRUE")
  B <- fit$params$B
  r = ncol(B)
  N = nrow(X)
  
  tmp.Estep <- E_step_su(X, Y, Z, fit$par, b_on_h = b_on_h)

  tmp.Estep <- within(tmp.Estep,{
                    mu_H = mu_U - mu_T %*% B
                    Mtt = ET[1:r,1:r,drop=FALSE]
                    Muu = ET[(r+1):(2*r), (r+1):(2*r),drop=FALSE]
                    Mtu = ET[1:r, (r+1):(2*r),drop=FALSE]
                    })

  # (t,h)
  mu_TH <- with(tmp.Estep, cbind(mu_T, mu_H))
  
  Mhh <- with(tmp.Estep, Muu - t(Mtu)%*%B - t(B)%*%Mtu + t(B)%*%Mtt%*%B)
  Mth <- with(tmp.Estep, Mtu - Mtt%*%B)
  Mtthh <- with(tmp.Estep, rbind(cbind(Mtt, Mth),
                 cbind(t(Mth), Mhh)))
  
  Sth <- with(tmp.Estep, Stu - Stt%*%B)
  Stthh <- with(tmp.Estep, rbind(cbind(Stt, Sth),
                 cbind(t(Sth), Shh)))
  
  # ## 
  # # test
  # Stthh_test <- matrix(0, 2*r, 2*r)
  # for (i in 1:N){
  #   Stthh_test <- Stthh_test + Mtthh + crossprod(mu_TH[i, , drop=F])
  # }
  # all.equal(Stthh_test/N, Stthh)
  
  aabb <- with(fit$params, cbind(a,b))
  
  chunk_1 <- chunk_2 <- chunk_3 <- chunk_4 <- term_3 <- matrix(0, 2*r, 2*r)
  
  for(i in 1:N){
    mu_TH_i <- mu_TH[i, ,drop=F]
    Stthh_i <- Mtthh + crossprod(mu_TH_i)
    chunk_1 <- chunk_1 + Z[i,]^2 * Stthh_i
    chunk_2 <- chunk_2 + Z[i,]*(t(mu_TH_i)%*%aabb%*%Stthh_i + t(t(mu_TH_i)%*%aabb%*%Stthh_i) +
                                  drop(aabb%*%t(mu_TH_i))*(Mtthh - crossprod(mu_TH_i)))
    chunk_4 <- chunk_4 + 2*Stthh_i%*%t(aabb)%*%aabb%*%Stthh_i + 
      drop(mu_TH_i%*%t(aabb)%*%aabb%*%t(mu_TH_i))*(Mtthh - crossprod(mu_TH_i)) +
      sum(diag(crossprod(aabb)%*%Mtthh))*Stthh_i
    
    # for the third term
    if(i!=1){
      for(j in 1:(i-1)){
        mu_TH_j <- mu_TH[j, ,drop=F]
        Stthh_j <- Mtthh + crossprod(mu_TH_j)
        Esi <- Z[i,]*t(mu_TH_i) - Stthh_i%*%t(aabb)
        Esj <- Z[j,]*t(mu_TH_j) - Stthh_j%*%t(aabb)
        term_3 <- term_3 + Esi%*%t(Esj)
      }
    }
  }
  chunk_3 <- t(chunk_2)
  chunk <- chunk_1 - chunk_2 - chunk_3 + chunk_4
  
  # (N/fit$par$sig2G * Stthh) %>% solve %>% diag %>% raise_to_power(0.5)
  
  info_ab = (N/fit$par$sig2G * Stthh) - (1/fit$par$sig2G^2 * chunk) - (1/fit$par$sig2G^2 * (term_3+t(term_3)))
  sd_ab <- info_ab %>% solve %>% diag %>% abs %>% raise_to_power(0.5)
  return(list(a = sd_ab[1:r], b = sd_ab[r+1:r]))
}

#' @export
generate_params_bi <- function(X=NULL, Y=NULL, Z=NULL, p=NULL, q=NULL, r, rx, ry, alpha_x = 0.1, alpha_y = 0.1, alpha_tu = 0.1, 
                               B=1, a0=0,a=t(2),b=t(1), type=c('o2m','random','specify')){
  type=match.arg(type)
  if(type=="o2m"){
    p = ncol(X)
    q = ncol(Y)
    fit <- o2m(X, Y, r, rx, ry, stripped=TRUE)
    x_tp <- with(fit, cbind(Tt,U-Tt))
    ab <- as.vector(glm(as.factor(Z)~x_tp, family = 'binomial')$coef)
    # print(ab)
    return(with(fit, {
      list(
        W = W.,
        Wo = suppressWarnings(orth(P_Yosc.)),
        C = C.,
        Co = suppressWarnings(orth(P_Xosc.)),
        B = abs(cov(Tt,U)%*%MASS::ginv(cov(Tt)))*diag(1,r),
        SigT = cov(Tt)*diag(1,r),
        SigTo = sign(rx)*cov(T_Yosc)*diag(1,max(1,rx)),
        SigUo = sign(ry)*cov(U_Xosc)*diag(1,max(1,ry)),
        SigH = cov(H_UT)*diag(1,r),
        sig2E = (ssq(X)-ssq(Tt)-ssq(T_Yosc))/prod(dim(X)) + 0.01,
        sig2F = (ssq(Y)-ssq(U)-ssq(U_Xosc))/prod(dim(Y)) + 0.01,
        a0 = ab[1],
        a = t(ab[1+1:r]),
        b = t(ab[(1+r)+1:r])
      )}))
  }
  
  if(type=="random"){
    outp <- list(
      W = orth(matrix(rnorm(p*r), p, r)+1),
      Wo = suppressWarnings(sign(rx)*orth(matrix(rnorm(p*max(1,rx)), p, max(1,rx))+seq(-p/2,p/2,length.out = p))),
      C = orth(matrix(rnorm(q*r), q, r)+1),
      Co = suppressWarnings(sign(ry)*orth(matrix(rnorm(q*max(1,rx)), q, max(1,ry))+seq(-q/2,q/2,length.out = q))),
      B = diag(1, r),
      # B = diag(sort(runif(r,1,1),decreasing = TRUE),r), # set to 1 for simulation
      SigT = diag(sort(runif(r,1,3),decreasing = TRUE),r),
      SigTo = sign(rx)*diag(sort(runif(max(1,rx),1,3),decreasing = TRUE),max(1,rx)),
      SigUo = sign(ry)*diag(sort(runif(max(1,ry),1,3),decreasing = TRUE),max(1,ry)),
      a0=a0,
      a = a,
      b = b
    )
    outp$SigH = diag(alpha_tu/(1-alpha_tu)*(diag(outp$SigT%*%outp$B^2)), r) #cov(H_UT)*diag(1,r),
    outp$SigU <- with(outp, SigT %*% B^2 + SigH)
    return(with(outp, {
      c(outp,
        sig2E = alpha_x/(1-alpha_x)*(mean(diag(SigT)) + mean(diag(SigTo)))/p,
        sig2F = alpha_y/(1-alpha_y)*(mean(diag(SigT%*%B^2 + SigH)) + mean(diag(SigUo)))/q)
    }))
  }
  
  if(type=="specify"){
    set.seed(0)
    outp <- list(
      W = orth(matrix(rnorm(p*r), p, r)+1),
      Wo = suppressWarnings(sign(rx)*orth(matrix(rnorm(p*max(1,rx)), p, max(1,rx))+seq(-p/2,p/2,length.out = p))),
      C = orth(matrix(rnorm(q*r), q, r)+1),
      Co = suppressWarnings(sign(ry)*orth(matrix(rnorm(q*max(1,rx)), q, max(1,ry))+seq(-q/2,q/2,length.out = q))),
      B = diag(B, r),
      # B = diag(sort(runif(r,1,1),decreasing = TRUE),r), # set to 1 for simulation
      SigT = diag(3,r),
      SigTo = sign(rx)*diag(3,max(1,rx)),
      SigUo = sign(ry)*diag(3,max(1,ry)),
      a0=a0,
      a = a,
      b = b
    )
    outp$SigH = diag(alpha_tu/(1-alpha_tu)*(diag(outp$SigT%*%outp$B^2)), r) #cov(H_UT)*diag(1,r),
    outp$SigU <- with(outp, SigT %*% B^2 + SigH)
    return(with(outp, {
      c(outp,
        sig2E = alpha_x/(1-alpha_x)*(mean(diag(SigT)) + mean(diag(SigTo)))/p,
        sig2F = alpha_y/(1-alpha_y)*(mean(diag(SigT%*%B^2 + SigH)) + mean(diag(SigUo)))/q)
    }))
  }
}

#' @export
generate_data_bi <- function(N, params, distr = rnorm){
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
  a0= params$a0
  a = params$a
  b = params$b
  
  p = nrow(W)
  q = nrow(C)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  Gamma = rbind(cbind(W, matrix(0,p,r), Wo, matrix(0,p,ry)),
                cbind(matrix(0,q,r), C, matrix(0,q,rx), Co))
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
  
  EF <- cbind(scale(matrix(rnorm(N*p), N))*sqrt(sig2E), 
              scale(matrix(rnorm(N*q), N))*sqrt(sig2F))
  
  mu_z <- a0 + M[,1:2*r]%*%t(cbind(a-b%*%B,b))
  prob_z <- 1/(1+exp(-mu_z))
  # print(cbind(M[,1:2*r],prob_z))
  Z <- rbinom(N, 1, prob_z)
  return(list(dat = cbind(M %*% t(Gamma) + EF, Z), trueTU = M[,1:2*r], mu_z = mu_z))
}


#' @export
E_step_bi <- function(X, Y, Z, params, level = level, Nr.core=1){
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
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  
  # Numerical integration of the common part and store the results (all sample combined)
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
    logl = logl
  )
}

## function for numerical integration of Q_ab
GH_Q_ab <- function(beta, l_n=list_nodes_th,Z.=Z,l_com=list_com_sub,de=dnmt){
  N = length(Z.)
  # log probability of z=1 and z=0
  z1 <- lapply(l_n, function(e){
    a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
    m <- max(0, a_tmp)
    -(m + log(exp(-m) + exp(a_tmp-m)))
  })
  z0 <- lapply(l_n, function(e){
    a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
    m <- max(0, a_tmp)
    -(m + log(exp(-m) + exp(a_tmp-m)))
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
  Q_ab <- lapply(1:N, function(e) Reduce("+", Map("*", int_diff_z[[e]], l_com[[e]])))
  # Q_ab <- rapply(Q_ab, f=function(x) ifelse(is.nan(x),0,x), how="replace")
  Q_ab <- Map("/", Q_ab, de)
  # Add N samples
  Q_ab <- Reduce("+", Q_ab)
  return(Q_ab)
}

## function for numerical integration of gradient of Q_ab
GH_grd_ab <- function(beta, l_n=list_nodes_th,Z.=Z,l_com=list_com_sub,de=dnmt){
  N = length(Z.)
  # Q'_ab
  z1 <- lapply(l_n, function(e){
    a_tmp <- as.numeric(-(cbind(1,e) %*% t(beta)))
    m <- max(0, a_tmp)
    exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)) * cbind(1,e)
  })
  z0 <- lapply(l_n, function(e){
    a_tmp <- as.numeric(cbind(1,e) %*% t(beta))
    m <- max(0, a_tmp)
    - exp(a_tmp-m)/(exp(-m) + exp(a_tmp-m)) * cbind(1,e)
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
  # grd_ab <- rapply(grd_ab, f=function(x) ifelse(is.nan(x),0,x), how="replace")
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
    if(orth_x){
      dcmp_varxo <- svd_orthpart(X, params$W, params$sig2E, mu_T, Stt, rx)
      params$Wo = dcmp_varxo$V
      params$SigTo = diag(x=dcmp_varxo$D, nrow=rx)
    }
    if(orth_y){
      dcmp_varyo <- svd_orthpart(Y, params$C, params$sig2F, mu_U, Suu, ry)
      params$Co = dcmp_varyo$V
      params$SigUo = diag(x=dcmp_varyo$D, nrow=ry)
    }
    
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
#' @param Z vector of 1 and 0. 1 for case, 0 for control.
Su_PO2PLS_bi <- function(X, Y, Z, r, rx, ry, steps = 1e5, tol = 1e-6, level=level, Nr.core =1, init_param= c('random','o2m'),
                         orth_type = "SVD", random_restart = FALSE){
  # if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F") %in% names(init_param))) {cat('using old fit \n'); params <- init_param}
  # else {params <- generate_params(X, Y, Z, r, rx, ry, type = init_param)}
  
  # debug_list <- list()
  # z_inf <- c()
  
  # Check Z
  
  
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
  mu_z <- with(params, a0 + X%*%W%*%(t(a)-B%*%t(b)) + Y%*%C%*%t(b))
  z_p <- 1/(1+exp(-mu_z))
  pre_z <- ifelse(z_p > 0.5, 1,0)
  levs <- c(0,1)
  zfit <- table(factor(pre_z, levs), factor(Z, levs))
  acc <- caret::confusionMatrix(zfit)$overall['Accuracy']
  
  message("Nr steps was ", i)
  message("Negative increments: ", any(diff(logl[1:i+1]) < 0), "\n",
          "Smallest increments: ", min(diff(logl[1:i+1])),
          "; Last increment: ", signif(logl[i+1]-logl[i],4))
  message("Log-likelihood: ", logl[i+1])
  message("Accuracy of z fit: ", acc)
  outputt <- list(params = params, logl = logl[0:i+1][-1], mu_z=mu_z, zfit = zfit, acc=acc) #, debug = debug_list)
  class(outputt) <- "PO2PLS"
  return(outputt)
}
