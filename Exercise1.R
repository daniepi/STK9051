#### STK 9051 - Project exam I ####
## Run Exploratory.R first
## c.) 
dat=as.numeric(data)

### By greedy (local) search, try to optimize likelihood. Here we want to optimize the assignmnet of the states
### to every observation. Thus the neighborhood is defined as a vector of length n, where inputs specify state
### for given observation.
### Function to generate intial parameter values
generate_initial_par = function(K=K, mu, Sigma, prob){
  ### K - number of hidden states
  ### mu - vector of length K, specifying initial guess for normal mean
  ### Sigma - vector of length K, specifying initial guess for normal std
  ### prob - probability of staying at given state
  if(length(Sigma)!=K){stop('Length of Sigma is not valid!')}
  if(length(mu)!=K){stop('Length of mu is not valid!')}
  if(length(prob)!=1){warning('The first value of prob is taken!')}
  P = matrix((1-prob[1])/(K-1), K, K); diag(P) = prob[1]
  list(pi=rep(1/K,K), p = P, mu=mu, sigma=Sigma)
}

likHoodState = function(dat, param, states){
  n=length(dat)
  for(t in 1:n){
    if(t == 1){
      forward = 1/K*dnorm(dat[t], param$mu[states[t]], param$sigma[states[t]])
      log_forward = log(forward)
    }
    else{
      forward = param$p[states[t-1], states[t]]*dnorm(dat[t], param$mu[states[t]],
                                                      param$sigma[states[t]])
      log_forward = log_forward + log(forward)
    }
  }
  log_forward
}


### By greedy (local) search, try to optimize likelihood. Here we want to optimize the assignmnet of the states
### to every observation. Thus the neighborhood is defined as a vector of length n, where inputs specify state
### for given observation.
param=generate_initial_par(K=3,c(0,0,0),c(0.1,0.55,1),0.9)

greedy_stock = function(dat, param, maxit=100, tol=10**-8,seed=123){
  n = length(dat)
  ## Generate initial state assignment
  set.seed(seed)
  initial_state = sample(1:3, n, replace = T)
  #initial_state = as.numeric(states)
  ## Try 'smart' initial assignment
  #initial_state = apply(Gauss, 1, which.max)
  ## Initialize log-likelihood 
  ll = likHoodState(dat, param, initial_state)
  tot_likhood = ll
  state=initial_state
  ### Generate matrix for all data points, where each column represent Gaussian density under each state
  Gauss = mapply(FUN=function(mu, sigma){dnorm(dat, mu, sigma)}, param$mu, param$sigma)
  
  for(iter in 1:maxit){
    #state=initial_state
    for(i in 1:n){
      ll.tmp = ll
      old.state=state[i]
      new.state = sample((1:3)[-state[i]],1)
      if(i==1){add = (log(1/K*Gauss[i,new.state])-log(1/K*Gauss[i,old.state])) + log(param$p[new.state, state[i+1]])-
        log(param$p[old.state, state[i+1]])}
      else if(i < n){
        add = log(param$p[new.state,state[i+1]]*Gauss[i,new.state]) -
          log(param$p[old.state,state[i+1]]*Gauss[i,old.state]) + log(param$p[state[i-1],new.state]) -
          log(param$p[state[i-1],old.state])
      }
      else{
        add = log(Gauss[i,new.state])-log(Gauss[i,old.state]) + log(param$p[state[i-1],new.state]) -
          log(param$p[state[i-1],old.state])
      }
      if(add > 0){
        ll = ll.tmp + add
        state[i] = new.state
      }
    }
    if(abs(ll-tot_likhood[iter])<tol){
      return(list(states=state, lHood=tot_likhood))
    }
    tot_likhood[iter+1] = ll
  }
  return(list(states=state, lHood=tot_likhood))
}
greedyState_Rnd = greedy_stock(dat,param,1000)

### Plot
greedyState_Rnd$states = xts(greedyState_Rnd$states, order.by = index(data))
ind=index(data['1996'])
state=factor(greedyState_Rnd$states[ind], labels=c("Low vol.",'Normal mode','High vol.'))
x=data.frame(x=ind, y=data[ind], state=state)
#tikz('statesGreedy96.tex', height=6, width=8)
ggplot(data=x, aes(x=ind, y=Volume))+
  geom_line()+geom_point(aes(col=state))+
  ylab(expression(Delta(log(V[t]))))+xlab('')+
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y")
#dev.off()
### Compare state assignmnet
### Similar assignment as EM
#table(states, greedyState_Rnd$states)

## d.) SANN
sann_stock = function(dat, param, maxit=1000, tol=10**-4,seed=123, tau=100){
  n = length(dat)
  ## Generate initial state assignment
  set.seed(seed)
  initial_state = sample(1:3, n, replace = T)
  #initial_state = as.numeric(states)
  ## Try 'smart' initial assignment
  #initial_state = apply(Gauss, 1, which.max)
  ## Initialize log-likelihood 
  ll = likHoodState(dat, param, initial_state)
  tot_likhood = ll
  state=initial_state
  ### Generate matrix for all data points, where each column represent Gaussian density under each state
  Gauss = mapply(FUN=function(mu, sigma){dnorm(dat, mu, sigma)}, param$mu, param$sigma)
  for(iter in 1:maxit){
    #state=initial_state
    u=runif(n)
    for(i in 1:n){
      ll.tmp = ll
      old.state=state[i]
      new.state = sample((1:3)[-state[i]],1)
      if(i==1){add = (log(1/K*Gauss[i,new.state])-log(1/K*Gauss[i,old.state])) + log(param$p[new.state, state[i+1]])-
        log(param$p[old.state, state[i+1]])}
      else if(i < n){
        add = log(param$p[new.state,state[i+1]]*Gauss[i,new.state]) -
          log(param$p[old.state,state[i+1]]*Gauss[i,old.state]) + log(param$p[state[i-1],new.state]) -
          log(param$p[state[i-1],old.state])
      }
      else{
        add = log(Gauss[i,new.state])-log(Gauss[i,old.state]) + log(param$p[state[i-1],new.state]) -
          log(param$p[state[i-1],old.state])
      }
      if(add > 0){## If positive change, we adopt new solution
        ll = ll.tmp + add
        state[i] = new.state
      }else{## Else
        prob=exp(iter*add/tau)
        if(u[i]<prob){
          ll = ll.tmp + add
          state[i] = new.state
        }
      }
    }
    if(abs(ll-tot_likhood[iter])<tol){
      return(list(states=state, lHood=tot_likhood))
    }
    tot_likhood[iter+1] = ll
  }
  return(list(states=state, lHood=tot_likhood))
}
sannState_Rnd = sann_stock(dat,param,10000)
##
### Plot
sannState_Rnd$states = xts(sannState_Rnd$states, order.by = index(data))
ind=index(data['1996'])
state=factor(sannState_Rnd$states[ind], labels=c("Low vol.",'Normal mode','High vol.'))
x=data.frame(x=ind, y=data[ind], state=state)
#tikz('statesSann96.tex', height=6, width=8)
ggplot(data=x, aes(x=ind, y=Volume))+
  geom_line()+geom_point(aes(col=state))+
  ylab(expression(Delta(log(V[t]))))+xlab('')+
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y")
#dev.off()
### Compare state assignmnet
### Similar assignment as EM
#table(states, sannState_Rnd$states)

### Show that local search is trapped in local optimum.
iter=c(length(sannState_Rnd$lHood),length(greedyState_Rnd$lHood))
method=factor(rep(c(1,2),iter), labels = c('SANN','Local search'))
x=data.frame(x=c(1:iter[1],1:iter[2]), y=c(sannState_Rnd$lHood, greedyState_Rnd$lHood), method=method)
#tikz('trap.tex', height=6, width=8)
ggplot(data=x, aes(x=x, y=y, col=method))+
  geom_line(lwd=2)+
  ylab(expression(P(C,x,theta)))+xlab('Iter')
#dev.off()