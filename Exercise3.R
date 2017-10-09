### STK 9051 - Project exam I
## Run Exploratory.R, Exercise1.R, Exercise2.R first

### Try to implement likelihood such that we can optimize it directly through numerical optimizer.
### We can rewrite the forward/backward procedure into the following two steps using matrix notation

### 1.) Forward: f(t) = pi_1%*%G(K, x_1) for t= 1, where pi_1 is initial state distribution (P(C1=k)) and
### G(K) is KxK diagonal matrix where entries are Gaussian densities.
### For t=2,,,n : f(t) = f(t-1)%*%P%*%G(K, x_t). This step corresponds to both prediction and updating at the same
### time. 

## Compute log-likelihood
lHood_optim = function(dat, param, rescale){
  ### param is a vector of transformed parameters of length: 2xK +(K-1)x(K+1)
  ### Determine K
  K = -1+sqrt(2+length(param))
  n=length(dat)
  ##### Transform parameters to original form
  pi=exp(param[1:(K-1)])/(1+exp(param[1:(K-1)])) ## Initial state distribution
  pi[K] = 1-sum(pi)
  P = exp(param[K:(K**2-1)])/(1+exp(param[K:(K**2-1)])) ## Transition matrix
  P = matrix(P,ncol=K-1)
  P = cbind(P, 1-rowSums(P))
  mu = param[K**2:(K**2+K-1)] ## mu
  if(rescale){mu=mu*mean(dat)}
  sigma=exp(tail(param, K)) ## sigma
  #####
  # In order to make it run faster, we only do forward step
  #####
  for(t in 1:n){
    if(t == 1){
      forward = pi%*%Gauss_mat(dat[t], mu=mu, sigma=sigma)
      sum_forward = sum(forward)
      log_forward = log(sum_forward)
      forward=forward/sum_forward
    }
    else{ # rescaling to avoid overfloat
      forward = forward%*%P%*%Gauss_mat(dat[t], mu=mu, sigma=sigma)
      sum_forward = sum(forward)
      log_forward = log_forward + log(sum_forward)
      forward=forward/sum_forward
    }
  }
  log_forward ## log-likelihood
}

func_Optimize = function(dat, param, rescale=F, method='BFGS', maxit=20){
  ## Here params are at the original form in a list i.e. :
  ## param$mu - vector of length K with Gaussian means
  ## param$sigma - vector of length K with Gaussian standard deviations
  ## param$p - KxK transition matrix
  ## param$pi - vector of length K with initial state distribution
  K=length(param$pi)
  pi = param$pi; P=param$p; mu=param$mu; sigma=param$sigma
  ### In order to go from [0,1] --> (-Inf,Inf) we do logit transform and remove the last value (overparametrized).
  pi_opt = log(pi/(1-pi))[-K] ## Here we might get problem with pi=1
  ### The same for transition matrix. We remove last column (overparametrized).
  P_opt = as.vector(log(P/(1-P))[,-K]) ## Here we might get problem with pi=1
  ### We take log transform of sigma: (0, Inf) --> (-Inf, Inf)
  sigma_opt = log(sigma)
  ### We might rescale mu for faster convergence.
  if(rescale){scale=mean(dat); mu_opt = mu/scale}
  else{mu_opt=mu}
  param_opt = c(pi_opt, P_opt, mu_opt, sigma_opt)
  #### Now the optimization via optim. Do it as maximization
  opt = optim(param_opt, fn=lHood_optim, dat=dat, rescale=rescale, method=method, 
              control=list(maxit=maxit, trace=3, fnscale=-1))
  pars=opt$par; val=opt$value
  ##### Transform parameters to original form
  pi=exp(pars[1:(K-1)])/(1+exp(pars[1:(K-1)])) ## Initial state distribution
  pi[K] = 1-sum(pi)
  P = exp(pars[K:(K**2-1)])/(1+exp(pars[K:(K**2-1)])) ## Transition matrix
  P = matrix(P,ncol=K-1)
  P = cbind(P, 1-rowSums(P))
  mu = pars[K**2:(K**2+K-1)] ## mu
  if(rescale){mu=mu*mean(dat)}
  sigma=exp(tail(pars, K)) ## sigma
  return(list(pi=pi, p=P, mu=mu, sigma=sigma, lHood= val))
}


##### Example
initial_param = generate_initial_par(K=3,c(0,0,0),c(0.1,0.55,1),0.9)
a=func_Optimize(dat, initial_param, rescale = T)
q.smo = HMM.E(data,a)$q.smo
statesOptim=apply(q.smo, 1, which.max)
statesOptim=xts(statesOptim, order.by = index(data))
table(statesOptim)
tapply(data, statesOptim, stat_fun)
### Plot the states on top of the series on some chosen subset
ind=index(data['1996'])
x=data.frame(x=ind, y=data[ind], state=factor(statesOptim[ind], labels=c("Low vol.",'Normal mode','High vol.')))
#tikz('Optimstates96.tex', height=6, width=8)
ggplot(data=x, aes(x=ind, y=Volume))+
  geom_line()+geom_point(aes(col=state))+
  ylab(expression(Delta(log(V[t]))))+xlab('')+
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y")
#dev.off()

initial_param = generate_initial_par(K=2,c(0,0),c(0.1,0.7),0.9)
b=func_Optimize(dat, initial_param, rescale = T)
###
table(states, statesOptim)
