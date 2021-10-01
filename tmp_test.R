test_t <- seq(-6,6,length=30)
den <- c()
lden <- c()

for(i in 1:30){
  den[i] <- with(params, exp(-0.5*(
    p*log(2*pi) + p*log(sig2E) +
      1/sig2E*((X[1, ]-test_t[i]*t(W))%*%t(X[1, ]-test_t[i]*t(W)))))) %>% as.numeric()
}

for(i in 1:30){
  lden[i] <- with(params, -0.5*(
    p*log(2*pi) + p*log(sig2E) +
      1/sig2E*((X[1, ]-test_t[i]*t(W))%*%t(X[1, ]-test_t[i]*t(W))))) %>% as.numeric()
}



plot(test_t, den)
plot(test_t, lden)
