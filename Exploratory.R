#### STK9051 -- Project 1 ####
## Exercise 1 ##
library("tseries", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4.1")
library("xts", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4.1")
library("PerformanceAnalytics", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4.1")
library("ggplot2", lib.loc="/mn/sarpanitu/modules/packages/R/3.4.2-prerelease-gcc/lib64/R/library")
library("tikzDevice", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4.1")
setwd("~/PhD/Autumn 2017/STK9051/Project")
sp500 = as.xts(get.hist.quote('^gspc', start='1950-01-01', end='2012-09-09', quote = 'Volume', compression='d'))
##
## Take log of first differences
data = diff(log(sp500))[-1]

### a.) plot both volume={V_t} and log-volume change={r_t}
#tikz('spVolume.tex', height = 6, width=8)
plot(sp500/10**9, grid.col='grey90', yaxis.right=F, lwd=0.2, main=paste('SP500 volume series in', expression(1e9)),
     ylab=expression(V[t]), grid.ticks.on='years', major.ticks='years')
#dev.off()
#tikz('spLogVolume.tex', height = 6, width=8)
plot(data, grid.col='grey90', yaxis.right=F, lwd=0.5, main='SP500 volume log-change',
     ylab=expression(V[t]), grid.ticks.on='years', major.ticks='years')
#dev.off()
### Clearly, pure volume measure departs far from stationary times series. Both the mean and variance are not stationary.
### The series of log-change looks much closer to white noise. Mean looks stationary, however the variance seems to
### vary with a time. In addition, there are a few peaks and lows which deviate from typical Gaussian distribution.
skew1 = skewness(sp500); kurt1 = PerformanceAnalytics::kurtosis(sp500, method = 'excess')
skew2 = skewness(data); kurt2 = PerformanceAnalytics::kurtosis(data, method = 'excess')
#expr=as.expression(mapply(FUN=function(x,i){bquote(hat(gamma[.(i)]) == .(x))}, round(c(skew1, kurt1),2), c(1,2)))
expr = c(paste("gamma[1] == ", round(skew1,2)),paste("gamma[2] == ", round(kurt1,2)))
#tikz('densVolume.tex', height = 6, width=8)
ggplot(data=sp500, aes(x=sp500/10**9))+
  geom_line(stat='density')+
  labs(title='Kernel density of SP500 volume series')+
  xlab(expression(V[t]))+
  theme(axis.title=element_text(size=14,face="bold"))+
  annotate("text", x = 9, y = c(4,3.5), label = expr, parse=T, cex=5)+
  annotate("point", x = 8, y = c(4,3.5))
#dev.off()
#plot(density(sp500/10**9), xlab=expression(V[t]), main='Kernel density of SP500 volume series')
#legend('topright', legend=expr, pch=c(16,16),col=rep(1,2),bty='n', cex=1.5)
expr = c(paste("gamma[1] == ", round(skew2,2)),paste("gamma[2] == ", round(kurt2,2)))
#tikz('densLogVolume.tex', height = 6, width=8)
ggplot(data=data, aes(x=data))+
  geom_line(stat='density')+
  labs(title='Kernel density of SP500 volume log-change')+
  xlab(expression(Delta(log(V[t]))))+
  theme(axis.title=element_text(size=14,face="bold"))+
  annotate("text", x = 2, y = c(2,1.5), label = expr, parse=T, cex=5)+
  annotate("point", x = 1.5, y = c(2,1.5))
#dev.off()

## Even though density of log-change is symmetric, it is still heavy-tailed. None of the series is fit the Gaussian
## process. One would need to account for stochastic volatility.
## Check i.i.d assumption.
#tikz('acfVolume.tex', height = 6, width=8)
chart.ACF(sp500, maxlag = 20, main=expression(V[t]))
#dev.off()
#tikz('acfLogVolume.tex', height = 6, width=8)
chart.ACF(data, maxlag = 20, main=expression(v[t]))
#dev.off()
#tikz('acfLogVolume2.tex', height = 6, width=8)
chart.ACF(data**2, maxlag = 20, main=expression(v[t]))
#dev.off()
# Volatility clustering

## qq-plots
y <- quantile(sp500/10**9, c(0.25, 0.75))
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
ggplot(sp500, aes(sample = sp500/10**9)) + stat_qq() + geom_abline(slope = slope, intercept = int, col='red')+
  labs(title='QQ-plot for SP500 volume series')

y <- quantile(data, c(0.25, 0.75))
x <- qnorm(c(0.25, 0.75))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]
ggplot(data, aes(sample = data)) + stat_qq() + geom_abline(slope = slope, intercept = int, col='red')+
  labs(title='QQ-plot for SP500 volume log-change')