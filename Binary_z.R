library(magrittr)
r=1
sigma = matrix(c(4,1.9,1.9,3),2)
fun_h <- function(l){
  crossprod(l$n) * l$w
}

GH_Intl(fun_mu, dim=2*r, level=20, mu=rep(0,2*r), sigma=sigma, params = 3, plot_nodes=T)
