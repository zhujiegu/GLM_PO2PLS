library(magrittr)

N=500;p=10;q=10
r=3
rx=2
ry=2


x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
z <- rnorm(N)

params_true <- generate_params(x, y, z, r, rx, ry, alpha = 0.4, type='random')
params_true$a <- 1*t(3:1)
params_true$b <- 0*t(r:1)
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
  beta <- (with(fit_o2$params, glm(Z~X%*%W)) %>% coef)[-1]
  pre_zo2 <- with(fit_o2$params, X_t %*% W %*% beta)
  r2_o2[i] <- with(glm(Z_t~pre_zo2) %>% summary, 1 - deviance/null.deviance)
}

boxplot(r2_o2,r2_su, names = c("PO2PLS","supervised"))
boxplot(r2_su-r2_o2)
abline(h=0, col="red")

diff_100_50 <- r2_su-r2_o2



fit$debug %>% sapply(function(e) e$a[1]) %>% c(params_true$a[1],.) %>% plot
abline(h=params_true$a[1])
E_step(X, Y, Z, params_true) %>% str
fit$params %>% str
params_true %>% str
t(params_true$W) %*% fit$params$W





E_step(X, Y, Z, params_true) %>% str
PO2PLS::E_step(X, Y, params_true)[-(1:2)] %>% str
t(params_true$W) %*% fit_o2$params$W

beta <- (with(fit_o2$params, glm(Z~X%*%W)) %>% coef)[-1]
pre_zo2 <- with(fit_o2$params, X_t %*% W %*% beta)
with(glm(Z_t~pre_zo2) %>% summary, 1 - deviance/null.deviance)
