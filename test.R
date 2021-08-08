source('SuPO2PLS.R')
library(OmicsPLS)
library(magrittr)

#######
# set up
N=30;p=2;q=10
r=3
rx=1
ry=1

# generate data
x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
z <- rnorm(N)

alpha=0.1
params_true <- generate_params(x, y, z, r, rx, ry, alpha = alpha, type='random')
params_true$a <- t(c(0,0,1))
params_true$b <- t(c(0,0,1))
params_true$sig2G = with(params_true, as.numeric(alpha/(1-alpha)*(a %*% SigT %*% t(a) + b %*% SigU %*% t(b))))
# params_true$sig2G <- 0*0.001
# params_true$SigH <- diag(0.00001,r)
dat <- generate_data(N, params_true)

dim(dat)

X <- dat[,1:p]
Y <- dat[,(p+1):(p+q)]
Z <- dat[,p+q+1]

X %>% abs %>% mean
Z %>% abs %>% mean

# fit supervised PO2PLS
fit <- Su_PO2PLS(X, Y, Z, 2, rx, ry, steps = 1000, init_param = "o2m")


with(fit$params, X%*%W/sig2E) %>% abs %>% mean
with(fit$params, Z%*%a/sig2G) %>% abs %>% mean
with(fit$params, sig2G/sig2E)
with(params_true, sig2G/sig2E)

params_true %>% str
fit$params %>% str

t(params_true$W) %*% fit$params$W
pre_z <- with(fit$params, X%*%W%*%t(a) + Y%*%C%*%t(b))
# pre_z <- with(params_true, X%*%W%*%t(a) + Y%*%C%*%t(b))

cor(pre_z,Z)^2


# fit$debug %>% sapply(function(e) e$a[1]) %>% c(params_true$a[1],.) %>% plot
# abline(h=params_true$a[1])
# E_step(X, Y, Z, params_true) %>% str
# fit$params %>% str
# params_true %>% str
# 
# 
# pre_z <- with(fit$params, X%*%W%*%t(a) + Y%*%C%*%t(b))
# with(glm(Z~pre_z) %>% summary, 1 - deviance/null.deviance)
# 
fit_o2 <- PO2PLS::PO2PLS(X, Y, 2, rx, ry, steps = 1000, init_param = "o2m")
# 
# 
# PO2PLS::E_step(X, Y, params_true)[-1:2)] %>% str
t(params_true$W) %*% fit_o2$params$W
t(fit$params$W) %*% fit_o2$params$W
# 
# with(with(fit_o2$params, glm(Z~X%*%W) %>% summary), 1 - deviance/null.deviance)
