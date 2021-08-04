source('SuPO2PLS.R')
library(OmicsPLS)
library(magrittr)

#######
# set up
N=30;p=10;q=10
r=1
rx=1
ry=1

# generate data
x <- matrix(rnorm(p*N), N,p)
y <- matrix(rnorm(q*N), N,q)
z <- rnorm(N)

params_true <- generate_params_su(x, y, z, r, rx, ry, B=1, a=1,b=1,type='random')
params_true$a0 <- 0
params_true$alpha <- with(params_true, alpha <- cbind(a0,a,b))
dat <- generate_data_su(N, params_true)

dim(dat)

X <- dat[,1:p]
Y <- dat[,(p+1):(p+q)]
Z <- ifelse(dat[,p+q+1] > 0, 1,0)

params <- params_true  

GH_Intl(fun_S, dim=2*r, level=6, X[1,],Y[1,],Z[1], params)

