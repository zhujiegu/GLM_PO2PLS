# fast1 method fit 1 joint-comp glm-PO2PLS binary model on filtered x and y sequentially, 
# updating only one pair of regression coefficients 
Su_PO2PLS_bi_fast1 <- function(X, Y, Z, r, rx, ry, steps, tol, level, 
                              Nr.core =1, init_param= 'o2m', orth_type = 'SVD', 
                              random_restart = 'F'){
  if(r==1){
    return(Su_PO2PLS_bi(X, Y, Z, 1, rx, ry, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                        init_param= init_param, orth_type = orth_type, random_restart = random_restart))
  }
  if(r>1){
    print('More than 1 joint comp, fitting GLM-PO2PLS Gaussian model first')
    fit_tmp <- Su_PO2PLS(X, Y, Z, r, rx, ry, steps = 1000, tol = tol, init_param= init_param)
    fit_tmp$params <- fit_tmp$params[names(fit_tmp$params)!='sig2G']
    
    for(k in 1:r){
      # filter data
      x_k <- with(fit_tmp, X - latent_var$mu_T[,-k,drop=F] %*% t(params$W[,-k,drop=F]) -
                    latent_var$mu_To %*% t(params$Wo))
      y_k <- with(fit_tmp, Y - latent_var$mu_U[,-k,drop=F] %*% t(params$C[,-k,drop=F]) -
                    latent_var$mu_Uo %*% t(params$Co))
      print(paste('fitting binary model for component',k))
      if(all(c("W","Wo","C","Co","B","SigT","SigTo","SigUo","SigH","sig2E","sig2F","a0","a","b") %in% names(init_param))){
        message('using modified old fit as initial parameters \n')
        params_k <- init_param
        params_k$W <- params_k$W[,k,drop=F]
        params_k$C <- params_k$C[,k,drop=F]
        params_k$B <- params_k$B[k,k,drop=F]
        params_k$SigT <- params_k$SigT[k,k,drop=F] 
        params_k$SigH <- params_k$SigH[k,k,drop=F]
        params_k$SigU <- params_k$SigU[k,k,drop=F]
        params_k$a0 <- params_k$a0
        params_k$a <- params_k$a[,k,drop=F]
        params_k$b <- params_k$b[,k,drop=F]
        params_k$Wo <- params_k$Wo*0
        params_k$Co <- params_k$Co*0
        params_k$SigTo <- matrix(0,1,1)
        params_k$SigUo <- matrix(0,1,1)
        
        fit_k <- Su_PO2PLS_bi(x_k, y_k, Z, 1, 0, 0, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                              init_param= params_k, orth_type = orth_type, random_restart = random_restart)
      }else{
        fit_k <- Su_PO2PLS_bi(x_k, y_k, Z, 1, 0, 0, steps = steps, tol = tol, level=level, Nr.core =Nr.core, 
                              init_param= init_param, orth_type = orth_type, random_restart = random_restart)
      }
      
      # update initial estimate
      fit_tmp$params$W[,k] <- fit_k$params$W
      fit_tmp$params$C[,k] <- fit_k$params$C
      fit_tmp$params$B[k,k] <- fit_k$params$B %>% as.numeric()
      fit_tmp$params$SigT[k,k] <- fit_k$params$SigT %>% as.numeric()
      fit_tmp$params$SigH[k,k] <- fit_k$params$SigH %>% as.numeric()
      fit_tmp$params$SigU[k,k] <- fit_k$params$SigU %>% as.numeric()
      fit_tmp$params$a0 <- fit_k$params$a0
      fit_tmp$params$a[,k] <- fit_k$params$a
      fit_tmp$params$b[,k] <- fit_k$params$b
      fit_tmp$latent_var$mu_T[,k] <- fit_k$latent_var$mu_T
      fit_tmp$latent_var$mu_U[,k] <- fit_k$latent_var$mu_U
    }
    return(fit_tmp)
  }
}
