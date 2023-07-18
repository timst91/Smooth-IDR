n=100


y_test=seq(0,60,l=40)
nu=Inf


nsim=100

H_init=seq(0.5,5,l=5)
#H_init=c(0.01,0.5,1,2.5,5)

C=c(0.25,0.5,1,4)

res_logS=matrix(nrow=length(C),ncol=length(H_init))
dens=list()

col=2:6


labels <- sapply(H_init, function(h) substitute(paste("h"[init], " = ", s), list( s= h)))


par(mfrow=c(2,2))

for (i in 1:length(C)){
  c=C[i]
  
  
  if(i==1){
    plot(y_test,dgamma(y_test,shape = sqrt(x_test),scale=min(max(x_test,2),8)),
         type='l',main=paste("c=",c),ylab="",xlab="Threshold",
         ylim=c(0,0.1))
    legend("topright",legend=labels,
                  col=2:6,lty=1,cex=.6,bty="n")
  }else{plot(y_test,dgamma(y_test,shape = sqrt(x_test),scale=min(max(x_test,2),8)),
             type='l',main=paste("c=",c),ylab="",xlab="Threshold",
             ylim=c(0,0.07))}
  
  for(j in 1:length(H_init)){
    density=numeric(length(y_test))
    print(c(i,j))
    h_init=H_init[j]
    for(sim in 1:nsim){
      X=runif(n,0,10)
      y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))
      
      if(sim%%25==0){print(sim)}
      local=smooth_IDR_density_h_opt(y,X,x_test=x_test,c=c,
                                     h_init = h_init,
                                     y_test = y_test,nu=Inf,
                                     progress = TRUE,
                                     normalize = 1)
      density=density+1/nsim*local$density
      
    }
    dens=append(dens,list(density))
    lines(y_test,density,col=col[j],lwd=.4)
  }
}

par(mfrow=c(1,1))
