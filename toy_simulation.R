source('SuPO2PLS.R')
library(magrittr)

N=50;p=100;q=100
r=2
rx=1
ry=1


x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
z <- rnorm(N)

params_true <- generate_params(x, y, z, r, rx, ry, alpha = 0.5, type='random')
params_true$a <- 1*t(r:1)
params_true$b <- 1*t(r:1)
params_true$sig2G <- with(params_true, as.numeric(0.1/(1-0.1)*(a %*% SigT %*% t(a) + b %*% SigU %*% t(b))))
params_true %>% str

################
r2_su <- r2_o2 <- c()
for(i in 1:100){
  dat <- generate_data(N+500, params_true)
  
  dim(dat)
  
  X <- dat[1:N,1:p] %>% scale(scale = F)
  Y <- dat[1:N,(p+1):(p+q)] %>% scale(scale = F)
  Z <- dat[1:N,p+q+1] %>% scale(scale = F)
  X_t <- dat[-(1:N),1:p] %>% scale(scale = F)
  Y_t <- dat[-(1:N),(p+1):(p+q)] %>% scale(scale = F)
  Z_t <- dat[-(1:N),p+q+1] %>% scale(scale = F)
  
  
  fit <- Su_PO2PLS(X, Y, Z, r, rx, ry, steps = 3000, init_param = "o2m")
  pre_z <- with(fit$params, X_t%*%W%*%t(a) + Y_t%*%C%*%t(b))
  r2_su[i] <- with(glm(Z_t~pre_z) %>% summary, 1 - deviance/null.deviance)
  
  
  fit_o2 <- PO2PLS::PO2PLS(X, Y, r, rx, ry, steps = 3000, init_param = "o2m")
  beta <- (with(fit_o2$params, glm(Z~X%*%W + Y%*%C)) %>% coef)[-1]
  pre_zo2 <- with(fit_o2$params, cbind(X_t%*%W, Y_t%*%C) %*% beta)
  r2_o2[i] <- with(glm(Z_t~pre_zo2) %>% summary, 1 - deviance/null.deviance)
}

boxplot(r2_o2,r2_su, names = c("PO2PLS","supervised"))
boxplot(r2_su-r2_o2)
