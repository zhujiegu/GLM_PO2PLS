library(OmicsPLS)
library(magrittr)
library(tidyverse)
library(pryr)
library(glmnet)
source('glmPO2PLS.R')

N=100
p = 500
q = 10
alpha_x = 0.1
alpha_y = .1
alpha_tu = 0.4
r=1
rx=1
ry=1          
B = 1
a0=0
a=rep(1,r) %>% t
b=rep(2,r) %>% t

params_true <- generate_params_bi(NULL,NULL,NULL,p,q,r, rx, ry, alpha_x,alpha_y,
                                  alpha_tu,B,a0,a,b,type='specify')

dat <- generate_data_bi(N, params_true)
dat_t <- generate_data_bi(1000, params_true)
trueTU <- dat$trueTU
mu_z <- dat$mu_z
dat <- dat$dat
trueTU_t <- dat_t$trueTU
mu_z_t <- dat_t$mu_z
dat_t <- dat_t$dat
X <- dat[,1:p] %>% scale(scale = F)
Y <- dat[,(p+1):(p+q)] %>% scale(scale = F)
Z <- dat[,p+q+1]

X_t <- dat_t[,1:p] %>% scale(scale = F)
Y_t <- dat_t[,(p+1):(p+q)] %>% scale(scale = F)
Z_t <- dat_t[,p+q+1]

fit_su <- glm_PO2PLS(X, Y, Z, r, rx, ry, 'Gaussian', steps = 1000, tol = 1e-5, level=45, 
                     Nr.core =2, init_param = params_true)
fit_bi <- glm_PO2PLS(X, Y, Z, 1, 0, 0, 'binomial', steps = 10, tol = 0.1, level=5, Nr.core =4, init_param = 'o2m')

fit_fast1 <- Su_PO2PLS_bi_fast1(X, Y, Z, r, rx, ry, steps=10, tol=0.1, level=5, 
                               Nr.core =4, init_param = params_true, orth_type = 'SVD', 
                               random_restart = FALSE)

############################
## Comp order and sign switch
# joint part
od <- apply(abs(crossprod(fit$params$W, params_true$W)),2,which.max)
sn <- diag(ifelse(diag(crossprod(fit$params$W[, od], params_true$W))>0, 1,-1), nrow = r)

fit$params$W <- (fit$params$W[, od] %>% as.matrix(drop=F))%*%sn
fit$params$C <- (fit$params$C[, od] %>% as.matrix(drop=F))%*%sn
fit$params$B <- fit$params$B[od, od]
fit$params$SigT <- fit$params$SigT[od, od] %>% as.matrix(drop=F)
fit$params$SigH <- fit$params$SigH[od, od] %>% as.matrix(drop=F)
fit$params$SigU <- fit$params$SigU[od, od] %>% as.matrix(drop=F)
fit$params$a <- fit$params$a[, od, drop = F] %*% sn
fit$params$b <- fit$params$b[, od, drop = F] %*% sn

# specific part
od <- apply(abs(crossprod(fit$params$Wo, params_true$Wo)),2,which.max)
sn <- diag(ifelse(diag(crossprod(fit$params$Wo[, od], params_true$Wo))>0, 1,-1), nrow = rx)
fit$params$Wo <- (fit$params$Wo[, od] %>% as.matrix(drop=F))%*%sn

od <- apply(abs(crossprod(fit$params$Co, params_true$Co)),2,which.max)
sn <- diag(ifelse(diag(crossprod(fit$params$Co[, od], params_true$Co))>0, 1,-1), nrow = ry)
fit$params$Co <- (fit$params$Co[, od] %>% as.matrix(drop=F))%*%sn

###########################
# measures
crossprod(params_true$W, fit_su$params$W)
crossprod(params_true$C, fit_su$params$C)
crossprod(params_true$Wo, fit_su$params$Wo)
crossprod(params_true$Co, fit_su$params$Co)

crossprod(params_true$W, fit_bi$params$W)
crossprod(params_true$C, fit_bi$params$C)
crossprod(params_true$Wo, fit_bi$params$Wo)
crossprod(params_true$Co, fit_bi$params$Co)

crossprod(fit_su$params$W, fit_fast1$params$W)
crossprod(fit_su$params$C, fit_fast1$params$C)
crossprod(fit_su$params$Wo, fit_fast1$params$Wo)
crossprod(fit_su$params$Co, fit_fast1$params$Co)

crossprod(params_true$W, fit_fast1$params$W)
crossprod(params_true$C, fit_fast1$params$C)
crossprod(params_true$Wo, fit_fast1$params$Wo)
crossprod(params_true$Co, fit_fast1$params$Co)

fit_fast1$params$a0; fit_fast1$params$a; fit_fast1$params$b
