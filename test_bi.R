source('SuPO2PLS.R')
source('SuPO2PLS_bi.R')
source('dat_generator.R')
source('Numerical_int.R')
library(OmicsPLS)
library(magrittr)
library(parallel)

#######
# set up
N=50;p=1000;q=10
r=1
rx=0
ry=0

# generate data
x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
z <- rnorm(N)

params_true <- generate_params_bi(x, y, z, r, rx, ry, alpha_x = 0.1, alpha_y = 0.1,
                                  alpha_tu = 0.3, a0=0,a=t(2),b=t(1),type='random')
dat <- generate_data_bi(N, params_true)

dim(dat)

X <- dat[,1:p]
Y <- dat[,(p+1):(p+q)]
Z <- dat[,ncol(dat)]

params <- params_true

# GH_Intl(fun_mu, dim=2*r, level=6, X[1,],Y[1,],Z[1], params)
# GH_Intl(fun_S, dim=2*r, level=6, X[2,],Y[2,],Z[2], params)
# GH_Intl(fun_1, dim=2*r, level=6, X[2,],Y[2,],Z[2], params)
GH_com(Nr.cores=1, level=6, X,Y,Z, params, plot_nodes=T)%>% 
  lapply(unlist) %>% lapply(sum) %>% unlist %>% log %>% sum
###############################################
# test numerical estimation of the log likelihood
# ! Manually modify the output of the fun_com function (exclude z to match PO2PLS)
true_l <- PO2PLS::E_step(X,Y,params)$logl
test_l <- c()
levels <- c(2:20,25,30,50)
for(l in levels){
  test_l <- append(test_l, GH_com(Nr.cores=4, level=l, X,Y,Z, params, plot_nodes=T) %>% 
                     lapply(unlist) %>% lapply(sum) %>% unlist %>% log %>% sum)
}
plot(x=levels, y = test_l)
abline(h = true_l, col = 'red')
###############################################

# fit supervised PO2PLS
fit <- Su_PO2PLS_bi(X, Y, Z, r, rx, ry, steps = 5, level=9, Nr.core =4,init_param = "o2m")
fit2 <- Su_PO2PLS_bi(X, Y, Z, r, rx, ry, steps = 10, level=30, Nr.core =4,init_param = params)
fit2 <- Su_PO2PLS_bi(X, Y, Z, r, rx, ry, steps = 10, level=30, Nr.core =1,init_param = 'random')

fit_o2 <- PO2PLS::PO2PLS(X, Y, r, rx, ry, steps = 200, init_param = "random")

fit_su <- Su_PO2PLS(X, Y, Z, r, rx, ry, steps = 200, init_param = "random")


crossprod(fit$params$W, params$W)
crossprod(fit$params$Wo, params$Wo)
crossprod(fit$params$C, params$C)
crossprod(fit$params$Co, params$Co)

crossprod(fit2$params$W, params$W)
crossprod(fit2$params$Wo, params$Wo)
crossprod(fit2$params$C, params$C)
crossprod(fit2$params$Co, params$Co)


params %>% str
fit$params %>% str
fit$logl %>% plot
fit$logl %>% plot(ylim=c(-170,-140))

#############################################
