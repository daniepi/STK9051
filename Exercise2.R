### STK 9051 - Project Exam I ###
## Run Exploratory.R, Exercise1.R first
#Hidden Markov model function for class-probabilities

HMM.E = function(x,param)
{
  #Assumes 
  #x is a vector of length n giving observations
  #param is a list containing
  #pi: The class probabilities at time 1
  #p: A KxK matrix giving the transition probabilities, p[l,k]=p(k|l)
  #mu: A K-vector giving expectations within each class
  #sigma: A K-vector giving the standard deviations within each class
  K = length(param$pi)
  n = length(x)
  q.upd = matrix(nrow=n,ncol=K)
  q.pred = matrix(nrow=n,ncol=K)
  #Initialization
  q.upd[1,] = param$pi*dnorm(x[1],param$mu,param$sigma)
  q.upd[1,] = q.upd[1,]/sum(q.upd[1,])
  #Forward equations
  for(i in 2:n){
    #Prediction
    for(k in 1:K){
      q.pred[i,k] = sum(param$p[,k]*q.upd[i-1,])
    }
    q.pred[i,] = q.pred[i,]/sum(q.pred[i,])
    #Updating
    q.upd[i,] = q.pred[i,]*dnorm(x[i],param$mu,param$sigma)
    q.upd[i,] = q.upd[i,]/sum(q.upd[i,])
  }
  #Backward equations
  q.smo = matrix(nrow=n,ncol=K)
  q.smo[n,] = q.upd[n,]
  for(i in (n-1):1)
  {
    for(k in 1:K)
      q.smo[i,k] = sum(q.upd[i,k]*param$p[k,]*q.smo[i+1,]/q.pred[i+1,]) 
  }
  
  #Sequence prob
  q.seq = array(NA,c(n,K,K))
  for(l in 1:K){
    for(k in 1:K){
      q.seq[2:n,l,k] = q.smo[2:n,k]*q.upd[1:(n-1),l]*param$p[l,k]/q.pred[2:n,k]
    }
  }
  list(q.pred=q.pred,q.upd=q.upd,q.smo=q.smo,q.seq=q.seq)
}


M_step = function(dat, param){
  tmp = HMM.E(x=dat, param=param); n=length(dat); K = length(param$pi)
  ## A list with 3 N*K matrices (prediction, update, smoothing)
  ## Plus an array (K -- N*K matrices) with sequence probabilities
  pi1 = tmp$q.smo[1,] ## Initial class distribution at time 1
  ## Now K*K transition matrix P
  P_tmp = matrix(0,K,K)
  for(k in 1:K){
    for(l in 1:K){
      P_tmp[l,k] = sum(tmp$q.seq[2:n,l,k])/sum(tmp$q.smo[1:(n-1),l])
    }
  }
  ## Now mu and sigma vectors for Gaussian densities
  mu_tmp = apply(tmp$q.smo, 2, FUN=function(x){sum(x*as.numeric(dat))})/colSums(tmp$q.smo)
  sigma2_tmp = rep(0,K)
  for(i in 1:K){
    sigma2_tmp[i] = sum(tmp$q.smo[,i]*(dat-mu_tmp[i])**2)/sum(tmp$q.smo[,i])
  }
  list(pi=pi1, p=P_tmp, mu=mu_tmp, sigma = sqrt(sigma2_tmp), q=tmp)
}
Gauss_mat = function(x, mu, sigma){
  diag(dnorm(x, mu, sigma))
}


## 2c.)
lhood_func = function(dat, param, q.pred){
  n=length(dat); K=length(param$pi)
  ff=rep(0,n)
  ## Sequentially add values to likelihood
  for(t in 1:n){
    tmp_f = 0
    if(t == 1){ff[t]=ff[t]+log(sum(param$pi*dnorm(dat[t], param$mu, param$sigma)))}
    else{
      tmp_f = tmp_f + sum(q.pred[t,]%*%Gauss_mat(dat[t], param$mu, param$sigma))
      ff[t] = ff[t-1] + log(tmp_f)
    }
  }
  ff
}

EM_algo = function(dat, init_param, numit=10, tol=10**-3){
  dat=as.numeric(dat)
  param=init_param
  K=length(param$pi); n=length(dat)
  lHood = 0
  ### EM iterations
  for(iter in 1:numit){
    tmp_est = M_step(dat, param)
    ff=lhood_func(dat, param=tmp_est, q.pred = tmp_est$q$q.pred)
    lHood[iter] = ff[n]
    if(iter > 1){
      if(abs(lHood[iter]-lHood[iter-1]) <= tol){
        print('Converged')
        return(list(params=tmp_est, lHood=lHood, lHood_t = ff))
      }
    }
    param=tmp_est
  }
  return(list(params=tmp_est, lHood=lHood, lHood_t = ff))
}

### Function to calculate sample statistics
stat_fun = function(x){
  res=c(mean(x),sd(x),PerformanceAnalytics::skewness(x), PerformanceAnalytics::kurtosis(x, method='excess'))
  names(res) = c('Mean','Std','Skewness','Kurtosis')
  round(res,3)
}


################# Examples ##################
################################
### K=3
initial_param = generate_initial_par(K=3,c(0,0,0),c(0.1,0.55,1),0.9)
aa = EM_algo(data, initial_param, numit=50, tol=10**-3)
### Likelihood as a function of EM iterations
toplot = data.frame(lhood=aa$lHood, iter=1:length(aa$lHood))
#tikz('EMlikhood.tex', height=6, width=8)
ggplot(data=toplot, aes(x=iter, y=lhood))+geom_line()+
  ylab(expression(l(theta)))
#dev.off()
### Latent states and sample statistics
states=apply(aa$params$q$q.smo, 1, which.max)
states=xts(states, order.by = index(data))
table(states)
tapply(data, states, stat_fun)
### Plot the states on top of the series on some chosen subset
ind=index(data['1996'])
state=factor(states[ind], labels=c("Low vol.",'Normal mode','High vol.'))
x=data.frame(x=ind, y=data[ind], state=state)
#tikz('states96.tex', height=6, width=8)
ggplot(data=x, aes(x=ind, y=Volume))+
  geom_line()+geom_point(aes(col=state))+
  ylab(expression(Delta(log(V[t]))))+xlab('')+
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y")
#dev.off()

### Now try for K=2
initial_param = generate_initial_par(K=2,c(0,0),c(0.1,0.7),0.9)
### Now try the same with K=2
cc = EM_algo(data, initial_param, numit=50, tol=10**-1)
states2=apply(cc$params$q$q.smo, 1, which.max)
states2=xts(states2, order.by = index(data))
table(states2)
tapply(data, states2, stat_fun)
state=factor(states2[ind], labels=c("Normal mode",'High vol.'))
x=data.frame(x=ind, y=data[ind], state=state)
#tikz('2states96.tex', height=6, width=8)
ggplot(data=x, aes(x=ind, y=Volume))+
  geom_line()+geom_point(aes(col=state))+
  ylab(expression(Delta(log(V[t]))))+xlab('')+
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y")
#dev.off()

### How do we compare K=2 vs. K=3 models?
### First we need to determine how many parameters are being estimated.
### In transition matrix we need to determine Kx(K-1) probabilities. In initial distribution vector we need to compute
### (K-1) probabilities and 2xK parameters for Gaussian densities
### |param| = 2xK + (K-1) + Kx(K-1) = (K-1)x(K+1)+2xK
K=length(aa$params$pi)
np = (K-1)*(K+1)+2*K
AIC1 = -2*(tail(aa$lHood,1) - np)
K=length(cc$params$pi)
np = (K-1)*(K+1)+2*K
AIC2 = -2*(tail(cc$lHood,1) - np)
## According to AIC, we should prefer model with K=3.

### For the chosen 'best' model calculate likelihood in time(i) and plot it.
lhood_t = xts(aa$lHood_t, order.by = index(data))
## Full scale
par(mfrow=c(1,2))
#tikz('likHoodTime.tex', height=6, width=8)
chart.TimeSeries(data['1996/1997'], main=expression(Delta(log(V[t]))), ylab='', lwd=0.5,
                 period.areas = "1996-11-20::1996-12-05", period.color = 'darkolivegreen1')
Min=which.min(data['1996/1997']); Max=which.max(data['1996/1997'])
points(c(Min,Max),c(data['1996/1997'][c(Min,Max)]),pch=16,col='red',cex=1)
chart.TimeSeries(lhood_t['1996/1997'], main=expression(l(theta)), ylab='',
                 period.areas = "1996-11-20::1996-12-05", period.color = 'darkolivegreen1')
points(c(Min,Max),lhood_t['1996/1997'][c(Min,Max)],pch=16,col='red',cex=1)
#dev.off()
par(mfrow=c(1,1))
### Zoomed-in
par(mfrow=c(1,2))
#tikz('likHoodTime96.tex', height=6, width=8)
chart.TimeSeries(data['1996-11/1996'], main=expression(Delta(log(V[t]))), ylab='',
                 period.areas = "1996-11-27::1996-12-03", period.color = 'darkolivegreen1')
chart.TimeSeries(lhood_t['1996-11/1996'], main=expression(l(theta)), ylab='',
                 period.areas = "1996-11-27::1996-12-03", period.color = 'darkolivegreen1')
#dev.off()
par(mfrow=c(1,1))
## Here we can see that 'unexpected' points in a series (huge drops and huge increase) decrease likelihood, since
## these points 'do not fit' the global structure. The drop in likelihood is less drastic with more complicated model
## where K=3