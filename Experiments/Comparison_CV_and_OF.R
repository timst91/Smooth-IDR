n=c(seq(100,1400,l=60),seq(1500,3000,l=10))

nsim=10
time=c()
cv=c()
cv_time=c()
of=c()
of_time=c()
nu=2.5

pb=progress_bar$new(total=length(n)*nsim)
for (N in n){
  h=(log(N)/N)^(1/9)
  c=0
  c_time=0
  o=0
  o_time=0
  for (i in 1:nsim){
    y=c()
    X=sort(runif(N,0,10))
    for (x in X) {y=append(y,rgamma(1,shape=x,scale=1/sqrt(x)))}
    CV_int=CV_h(y,X,h,nu)
    OF_int=OF_h(y,X,h,nu)
    c=c+1/nsim*CV_int$logS
    c_time=c_time+1/nsim*CV_int$time
    o=o+1/nsim*OF_int$logS
    o_time=o_time+1/nsim*OF_int$time
    pb$tick()
  }
  print(c(o,c))
  cv=append(cv,c)
  cv_time=append(cv_time,c_time)
  of=append(of,o)
  of_time=append(of_time,o_time)
  
  
}
pb$terminate()


par(mfrow=c(1,2))
plot(n,cv,type='l',ylab="")
lines(n,of,col=2)
legend(x=1500,1.485, legend = expression(CV(h[n]), OF(h[n])), 
       col = c(1, 2), lty = 1, bty = "n", xpd = NA,
       x.intersp = 0.5)


plot(n,cv_time,type='l',ylim=c(min(of_time),max(cv_time)),ylab="",main = " Computational time (seconds)")
lines(n,of_time,col=2)
legend("topleft", legend = expression(CV(h[n]), OF(h[n])), 
       col = c(1, 2), lty = 1, bty = "n", xpd = NA,x.intersp = 0.5)

par(mfrow=c(1,1))
