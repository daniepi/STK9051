SMC_param = function(m=10**5, nT=length(Y), alpha, beta, sigma_y=1){
  ## Create objects
  tau.unique = rep(0, nT)
  tau_sim = matrix(0, nT, m)
  X = matrix(0, nT, m)
  N.eff = rep(0,nT)
  X.SIR.hat = X2.SIR.hat = rep(0,nT)
  CI = CI2 = matrix(0,nT,2)
  N.unique = matrix(0, nT, nT)
  var.weight=rep(0,nT)
  prob=c(0.025,0.975)
  ##########################
  ## First data point
  tau_sim[1,] = rgamma(m, alpha, beta) # prior
  X[1, ] = rnorm(m, 0, 10)
  #X[1, ] = rnorm(m, 0, 1/sqrt(tau_sim[1,]))
  w = dnorm(Y[1], X[1,]**2/25, sigma_y)
  w=w/sum(w)
  N.eff[1] = 1/sum(w**2)
  ## resample
  ind = sample(1:m, m, replace = T, prob = w)
  X[1,] = X[1,ind]
  tau_sim[1,] = tau_sim[1,ind]
  #w=rep(1/m,m)
  X.SIR.hat[1] = sum(X[1,])/m
  X2.SIR.hat[1] = sum(X[1,]**2)/m
  CI[1,] = quantile(X[1,], probs = prob)
  CI2[1,] = quantile(X[1,]**2, probs = prob)
  N.unique[1,1] = length(unique(X[1,]))
  tau.unique[1] = length(unique(tau_sim[1,]))
  var.weight[1] = var(w)
  
  mean_func = function(x, y){
    0.5*x+25*x/(1+x**2)+
      8*cos(1.2*(y-1))
  }
  mu = matrix(0,nT,m)
  add = rep(0, m)
  for(t in 2:nT){
    mu[t,] = mean_func(X[t-1,], t) # current mean
    #add = add + (X[t-1,]-mu[t-1,])**2 # update Gamma parameter
    #tau_sim[t,] = rgamma(m, alpha+t/2,  beta+0.5*add) # posterior
    X[t, ] = rnorm(m, mu[t,], sqrt(1/tau_sim[t-1,])) # state process
    add = add + (X[t,]-mu[t,])**2 # update Gamma parameter
    w = dnorm(Y[t], X[t,]**2/25, sigma_y) 
    w=w/sum(w)
    N.eff[t] = 1/sum(w**2)
    ## resample
    ind = sample(1:m, m, replace = T, prob = w)
    X[1:t,] = X[1:t,ind]
    add = add[ind]
    #tau_sim[t,] = tau_sim[t,ind]
    tau_sim[t,] = rgamma(m, alpha+(t-1)/2,  beta+0.5*add) # posterior
    #w=rep(1/m,m)
    X.SIR.hat[t] = sum(X[t,])/m
    X2.SIR.hat[t] = sum(X[t,]**2)/m
    CI[t,] = quantile(X[t,], probs = prob)
    CI2[t,] = quantile(X[t,]**2, probs = prob)
    N.unique[1:t,t] = apply(X[1:t,], 1, function(y){length(unique(y))})
    tau.unique[t] = length(unique(tau_sim[t,]))
    var.weight[t] = var(w)
  }
  N.unique[N.unique == 0] = NA
  res1 = data.frame(estimate = X.SIR.hat, CI = CI)
  res2 = data.frame(estimate = X2.SIR.hat, CI = CI2)
  return(list(X=X, var.weight=var.weight, N.unique=N.unique, N.eff=N.eff,
              res1=res1, res2=res2, tau=tau_sim, tau.unique = tau.unique))
}