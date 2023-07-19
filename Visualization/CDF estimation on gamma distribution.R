#install.packages("devtools")
library(devtools)
#loads the neccessary functions
source_url("https://raw.githubusercontent.com/timst91/Smooth-IDR/main/R/all.R")




######### Front page plot:
n=100


X=runif(n,0,10)
y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))

x_test=6

idr.fit=idr(y,data.frame(X))
pred.idr=predict(idr.fit,data=data.frame(X=x_test))
plot(pred.idr,lwd=.1)

#y_test=pred.idr[[1]]$points

y_test=seq(0,70,l=40)
nu=Inf

lines(y_test,smooth_IDR_CDF(y,X,nu=nu,x_test=x_test,y_test = y_test,
                            progress = TRUE),
      type='l',col=3)

lines(y_test,pgamma(y_test,shape = sqrt(x_test),scale=min(max(x_test,2),8)),
      col=2)

lines(y_test,smooth_IDR_CDF_h_opt(y,X,nu=nu,c=1,x_test=x_test,y_test = y_test,
                            progress = TRUE)$cdf,col=6)
lines(y_test,smooth_IDR_CDF_h_opt(y,X,nu=nu,c=0.5,x_test=x_test,y_test = y_test,
                                  progress = TRUE)$cdf,col=6,lty=2,lwd=.5)
lines(y_test,smooth_IDR_CDF_h_opt(y,X,nu=nu,c=2,x_test=x_test,y_test = y_test,
                                  progress = TRUE)$cdf,col=6,lty=3,lwd=.8)

legend("right",c(expression(F[x=6](y)),"IDR",
                       "smooth IDR","smooth IDR local bw c=0.5",
                       "smooth IDR local bw c=1", "smooth IDR local bw c=2"),
       lty=c(1,1,1,2,1,3),col=c(2,"blue",3,6,6,6),
       bty="n",text.width = 10)

##############################




col= c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 6)

labels=sapply(H_init, function(h) substitute(paste("h"[init], " = ", s), list( s= h)))
legend_items = matrix(labels[1:4], nrow = 2, byrow = TRUE)
legend_items=rbind(legend_items,c("",labels[5]))


par(mfrow=c(2,2))

for (i in 1:length(C)){
  c=C[i]
  
  
  if(i==1){
    plot(y_test,pgamma(y_test,shape = sqrt(x_test),scale=min(max(x_test,2),8)),
         type='l',main=paste("c=",c),ylab="",xlab="Threshold"
        )
   
    
  }else{plot(y_test,pgamma(y_test,shape = sqrt(x_test),scale=min(max(x_test,2),8)),
             type='l',main=paste("c=",c),ylab="",xlab="Threshold",
             )}
  
  for(j in 1:length(H_init)){
    cdf=numeric(length(y_test))
    print(c(i,j))
    h_init=H_init[j]
    for(sim in 1:nsim){
      X=runif(n,0,10)
      y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))
      
      if(sim%%25==0){print(sim)}
      local=smooth_IDR_CDF_h_opt(y,X,x_test=x_test,
                                     c=c,
                                     h_init = h_init,
                                     y_test = y_test,nu=Inf,
                                     progress = TRUE)
      cdf=cdf+1/nsim*local$cdf
      
    }
    cdfs=append(cdf,list(cdf))
    lines(y_test,cdf,col=col[j],lwd=.4)
  }
}



plot(y_test,pgamma(y_test,shape = sqrt(x_test),scale=min(max(x_test,2),8)),
     type='l',main="Global bandwidth CDF",ylab="",xlab="Threshold"
     )


labels2 =sapply(H_init, function(h) substitute(paste("h", " = ", s), list( s= h)))
legend_items = matrix(labels2[1:4], nrow = 2, byrow = TRUE)
legend_items=rbind(legend_items,c("",labels2[5]))


legend("bottomright",legend=legend_items,
       col=c(col[1],col[3],NA,col[2],col[4],col[5]),
       lty=1,cex=1,bty="n",ncol = 2)


for( i in (1:length(H_init))){
  density=numeric(length(y_test))
  for (sim in 1:nsim) {
    X=runif(n,0,10)
    y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))
    
    density=density+1/nsim*smooth_IDR_CDF(y,X,x_test = x_test,y_test = y_test,
                                              nu=Inf,h=H_init[i])
  }
  lines(y_test,density,lwd=.4,col=col[i])
}


par(mfrow=c(1,1))

