###########################################################
# Convergence Plots
####
# Plot a and b trends
a_plot1 <- sapply(fit$debug, function(e) e$aa[1,1])
a_plot2 <- sapply(fit$debug, function(e) e$aa[1,2])
b_plot1 <- sapply(fit$debug, function(e) e$bb[1,1])
b_plot2 <- sapply(fit$debug, function(e) e$bb[1,2])

start = 1; end = 1000
plot(a_plot1[start:end])
points(a_plot2[start:end], col = "red")
points(b_plot1[start:end], col="blue")
points(b_plot2[start:end], col="green")
# abline(h=params_true$a[start])

# Plot sig2G trends
sig2G_plot <- sapply(fit$debug, function(e) e$sig2G)
plot(sig2G_plot[start:end], col="purple")

# Plot determinant of Sttuu
cortu_plot <- sapply(fit$debug, function(e) det(e$Sttuu))
plot(cortu_plot[start:end], col="purple")

plot(sig2G_plot)
plot(fit$logl[start:end])

# Plot t(Z/N) %*% mu_T, mu_U
ZTU_plot <- sapply(fit$debug, function(e) e$ZTU)
plot(ZTU_plot[1,start:end])
plot(ZTU_plot[2,start:end], col = "red")
plot(ZTU_plot[3,start:end], col="blue")
plot(ZTU_plot[4,start:end], col="green")

  #############################################
  # browser inside Su_PO2PLS function
  
  E_next$EE %>% cor
  cor(E_next$mu_T, E_next$mu_U)
  cor(E_next$mu_U)
  cor(E_next$mu_T)
  
  E_next %>% str
  params_max %>% str
  params_max$W %>% cor
  
  
  with(fit$params, cor(X%*%W, Y%*%C))
  with(fit$params, cor(X%*%W))
  with(fit$params, cor(Y%*%C))
  
  E_next$Mtilde
  E_next$GammaEFG
  cor(E_next$GammaEFG %*% E_next$Mtilde)
  cor(E_next$EE)
  all.equal(cbind(X,Y,Z)%*%(E_next$GammaEFG %*% E_next$Mtilde),E_next$EE)
  cor(Y)
  
  with(E_next, (t(Z/N)%*%mu_U - params_next$a%*%t(Stu)) %*% MASS::ginv(Suu))
  with(E_next, params_next$a%*%t(Stu))
  
  
  MASS::ginv(E_next$Stt)
  
  params_next %>% str
  
  with(params_true, cor(X%*%W, Y%*%C))
  with(fit$params, cor(X%*%W, Y%*%C))
  