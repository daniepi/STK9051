MCMC = function(m=10**4, nT=length(Y), sigma_x=10, sigma_y=1, prop_sigma=1){
  ### Initialize X0 = 0
  X = matrix(0, m, nT)
  mean_func = function(x, y){
    0.5*x+25*x/(1+x**2)+
      8*cos(1.2*(y-1))
  }
  acc = rep(0,nT)
  prop_sigma_tmp = rep(prop_sigma,nT)
  
  ### MCMC step
  for(i in 2:m){
    if(i > m/10){
      w = acc==0
      if(any(w)){
        prop_sigma_tmp[!w] = prop_sigma
        prop_sigma_tmp[w] = prop_sigma_tmp[w]+5
        x.prop = X[i-1,] + rnorm(nT,0,prop_sigma_tmp)
      }
      else{
        x.prop = X[i-1,] + rnorm(nT,0,prop_sigma)
      }
    } else{
      x.prop = X[i-1,] + rnorm(nT,0,prop_sigma)
    }
    for(t in 1:nT){
      if(t==1){
        new = log(dnorm(x.prop[t],0,sigma_x)*dnorm(X[i-1,t+1],mean_func(x.prop[t], t+1), sigma_x)*dnorm(Y[t], x.prop[t]**2/25, sigma_y))
        old = log(dnorm(X[i-1,t],0,sigma_x)*dnorm(X[i-1,t+1],mean_func(X[i-1,t], t+1),sigma_x)*dnorm(Y[t], X[i-1,t]**2/25, sigma_y))
      } else if(t==nT){
        new = log(dnorm(x.prop[t],mean_func(X[i-1,t-1], t), sigma_x)*dnorm(Y[t], x.prop[t]**2/25, sigma_y))
        old = log(dnorm(X[i-1,t],mean_func(X[i-1,t-1], t), sigma_x)*dnorm(Y[t], X[i-1,t]**2/25, sigma_y))
      } else{
        new = log(dnorm(x.prop[t],mean_func(X[i-1,t-1], t), sigma_x)*dnorm(X[i-1,t+1],mean_func(x.prop[t], t+1), sigma_x)*dnorm(Y[t], x.prop[t]**2/25, sigma_y))
        old = log(dnorm(X[i-1,t],mean_func(X[i-1,t-1], t), sigma_x)*dnorm(X[i-1,t+1],mean_func(X[i-1,t], t+1), sigma_x)*dnorm(Y[t], X[i-1,t]**2/25, sigma_y))
      }
      if(old == -Inf & new != -Inf){
        X[i,t] = x.prop[t]
        acc[t] = acc[t]+1
        #cat("Update at",t,"\n")
      }
      else if(old == -Inf & new == -Inf){
        X[i,t] = X[i-1,t]
        #cat("Stuck at",t,'both',"\n")
      }
      else if(new == -Inf){
        X[i,t] = X[i-1,t]
        #cat("Stuck at",t,"\n")
      }
      else{
        if(runif(1) < exp(new-old)){
          X[i,t] = x.prop[t]
          acc[t] = acc[t]+1
        } else{
          X[i,t] = X[i-1,t]
        } 
      }
    }
  }
  return(list(X=X, acc=acc/m, w=which(acc==0)))
}