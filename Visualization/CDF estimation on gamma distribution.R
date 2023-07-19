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





y=c()
n=100
X=sort(runif(n,0,10))

for (x in X) {y=append(y,rgamma(1,shape=sqrt(x),scale=min(max(x,2),8)))}

x_test=5



idr_fit=idr(y,data.frame(X))

pred=predict(idr_fit,data=data.frame(X=x_test))

#y_test=pred[[1]]$points
y_test=seq(-1,37,l=500)

true_cdf=pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8))


NU=c(seq(2.001,50,l=100),Inf)
h=seq(0.1,15,l=50)

grid2=expand.grid(NU,h)

of_h_nu=c()
pb <- progress_bar$new(total=nrow(grid2))

for (i in 1:nrow(grid2)){
  
  nu=as.numeric(grid2[i,])[1]
  h=as.numeric(grid2[i,])[2]
  of=OF_h(y,X,h,nu)$logS
  of_h_nu=append(of_h_nu,
                 of)
  #print(of)
  pb$tick()
}

nu_opt_init=as.numeric(grid2[which.min(of_h_nu),])[1]
h_opt_init=as.numeric(grid2[which.min(of_h_nu),])[2]
pb$terminate()


c=seq(0.25,10,l=30)

grid=expand.grid(NU,c)

of_c_nu=c()
pb <- progress_bar$new(total=nrow(grid))

for (i in 1:nrow(grid)){
  
  nu=as.numeric(grid[i,])[1]
  c=as.numeric(grid[i,])[2]
  of_c_nu=append(of_c_nu,
                 OF_c(y,X,c,nu,nu_init=nu_opt_init,h_init=h_opt_init))
  pb$tick()
}

nu_opt=as.numeric(grid[which.min(of_c_nu),])[1]
c_opt=as.numeric(grid[which.min(of_c_nu),])[2]
pb$terminate()

smooth_CDF=smooth_IDR_CDF(y,X,h=h_opt_init,nu=nu_opt_init,y_test=y_test,x_test=x_test)
smooth_CDF_local=smooth_IDR_CDF_h_opt(y,X,h_init=h_opt_init,nu_init=nu_opt_init,
                                      c=c_opt,nu=nu_opt,
                                      y_test=y_test,x_test=x_test)
smooth_CDF_local_uncorrected=smooth_IDR_CDF_h_opt(y,X,h_init=h_opt_init,nu_init=nu_opt_init,
                                                  c=c_opt,nu=nu_opt,
                                                 y_test=y_test,x_test=x_test,
                                                  increasing = FALSE)

plot(pred)
lines(y_test,smooth_CDF,col=2)
lines(y_test,smooth_CDF_local,col=1)
lines(y_test,smooth_CDF_local_uncorrected,col=1,lty=2)
lines(y_test,true_cdf,lty=2,col="darkgreen")

Legend1= substitute(paste("smooth IDR with global bw, ", h[opt], " = ", H, " , ", nu[opt]," = ",nU)
                    , list(H = h_opt_init, nU=nu_opt_init))

Legend2=substitute(paste("smooth IDR with local bw, ", c[opt], " = ", C,
                         " , ", nu[opt]," = ",nU), list(C = c_opt, nU=nu_opt))



legend("right",c("IDR", Legend1, Legend2,
                 "smooth IDR with local bw, uncorrected for monotonicity",
                      expression(F[5](y))),
       col=c("blue",2,1,1,"darkgreen"),lty=c(1,1,1,2,2))






