SMC = function(m=10**5, nT=length(Y), sigma_x=10, sigma_y=1){
  ## Create objects
  X = XX = matrix(0, nT, m)
  N.eff = rep(0,nT)
  X.SIR.hat = X2.SIR.hat = rep(0,nT)
  CI = CI2 = matrix(0,nT,2)
  N.unique = matrix(0, nT, nT)
  var.weight=rep(0,nT)
  prob=c(0.025,0.975)
  ##########################
  ## First data point
  X[1, ] = XX[1,] = rnorm(m, 0, sigma_x)
  w = dnorm(Y[1], X[1,]**2/25, sigma_y)
  w=w/sum(w)
  var.weight[1] = var(w)
  N.eff[1] = 1/sum(w**2)
  ## resample
  ind = sample(1:m, m, replace = T, prob = w)
  X[1,] = X[1,ind]
  #w=rep(1/m,m)
  X.SIR.hat[1] = sum(X[1,])/m
  X2.SIR.hat[1] = sum(X[1,]**2)/m
  CI[1,] = quantile(X[1,], probs = prob)
  CI2[1,] = quantile(X[1,]**2, probs = prob)
  N.unique[1,1] = length(unique(X[1,]))
  mean_func = function(x, y){
    0.5*x+25*x/(1+x**2)+
      8*cos(1.2*(y-1))
  }
  for(t in 2:nT){
    X[t, ] = rnorm(m, mean_func(X[t-1,], t), sigma_x)
    w = dnorm(Y[t], X[t,]**2/25, sigma_y)
    w=w/sum(w)
    var.weight[t] = var(w)
    N.eff[t] = 1/sum(w**2)
    ## resample
    ind = sample(1:m, m, replace = T, prob = w)
    X[1:t,] = X[1:t,ind]
    XX[t,] = X[t,ind]
    #w=rep(1/m,m)
    X.SIR.hat[t] = sum(X[t,])/m
    X2.SIR.hat[t] = sum(X[t,]**2)/m
    CI[t,] = quantile(X[t,], probs = prob)
    CI2[t,] = quantile(X[t,]**2, probs = prob)
    N.unique[1:t,t] = apply(X[1:t,], 1, function(y){length(unique(y))})
  }
  N.unique[N.unique == 0] = NA
  res1 = data.frame(estimate = X.SIR.hat, CI = CI)
  res2 = data.frame(estimate = X2.SIR.hat, CI = CI2)
  return(list(X=X, var.weight=var.weight, N.unique=N.unique, N.eff=N.eff,
              res1=res1, res2=res2, XX=XX))
}