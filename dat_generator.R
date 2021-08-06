#' @export
generate_params_bi <- function(X, Y, Z, r, rx, ry, alpha_x = 0.1, alpha_y = 0.1,
                               alpha_tu = 0.1, a0=0,a=t(2),b=t(1), type=c('o2m','random')){
  type=match.arg(type)
  p = ifelse(is.matrix(X) | type != "random", ncol(X), X)
  q = ifelse(is.matrix(Y) | type != "random", ncol(Y), Y)
  if(type=="o2m"){
    stop("use random now")
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

  prob_z <- 1/(1+exp(-(a0 + M[,1:2*r]%*%t(cbind(a-b%*%B,b)))))
  # print(cbind(M[,1:2*r],prob_z))
  Z <- rbinom(N, 1, prob_z)
  cbind(M %*% t(Gamma) + EF, Z)
}