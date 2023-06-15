library(tictoc)
library(isodistrreg)

OF_h=function(y,x,h,nu){
  tic()
  sort_ind= match(sort(unique(x)), x)
  n=length(sort_ind)
  x=x[sort_ind]
  y=y[sort_ind]
  idr_fit=idr(y,data.frame(x))
  
  #pb <- progress_bar$new(total=n)
  if (nu!=Inf){
    K_h=function(u){1/h*(1/(nu*pi)^(0.5))*gamma((nu+1)/2)/gamma(nu/2)*(1+(u/h)^2/nu)^(-(nu+1)/2)}
  }
  else{
    K_h=function(u){1/h*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)}
  }
  
  cdf=idr_fit$cdf
  logS=0
  count=0
  wk=c()
  for (i in 1:n){
    
    w_1=cdf[i,c(1,1:(length(unique(y))-1))]
    w_1[1]=0
    w_i=cdf[i,]-w_1
    j_i=which(match(sort(unique(y)), y[i])==1)
    w_j_i=w_i[j_i]
    w_i[j_i]=0
    if(w_j_i!=1){
      w_i=w_i/(1-w_j_i)
    }else{
      count=count+1
      w_i=w_i}
    
    
    K=sapply(sort(unique(y)),function(u){K_h(u-y[i])})
    wik=w_i%*%K
    
    if(sum(wik)!=0){
      logS=logS-1/n*log(wik)
      wk=append(wk,wik)
      
    }else{
      count=count+1
      if(length(wk)>=1){
        wik=ifelse(length(wk)<=5,mean(wk),
                   mean(wk[(length(wk)-5):length(wk)]))
        
        logS=logS-1/n*log(wik)
        wk=append(wk,wik)
      }
      
      
    }
    
    #pb$tick()
    
  }
  
  time=toc(quiet = TRUE)
  #print(n/(n-count))
  return(list(logS=n/(n-count)*logS,time=unname(time$toc-time$tic)))
  
}


CV_h=function(y,x,h,nu,progress=FALSE){
  tic()
  sort_ind= match(sort(unique(x)), x)
  n=length(sort_ind)
  x_unique=x[sort_ind]
  logS=0
  
  if (progress){pb <- progress_bar$new(total=length(x_unique))}
  
  
  for (i in 1:length(x_unique)){
    model=smooth_IDR_density(y=y[sort_ind][-i],x=x_unique[-i],h=h,nu=nu,
                             y_test = y[sort_ind][i],x_test =x_unique[i])
    logS=logS-1/n*log(model)
    
    if(progress){pb$tick()}
  }
  # Close the progress bar
  time=toc(quiet=TRUE)
  
  return(list(logS=logS,time=unname(time$toc-time$tic)))
  
}

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
