y=c()
n=500
X=sort(runif(n,0,10))

for (x in X) {y=append(y,rgamma(1,shape=sqrt(x),scale=min(max(x,2),8)))}

x_test=5



idr_fit=idr(y,data.frame(X))

pred=predict(idr_fit,data=data.frame(X=x_test))

y_test=pred[[1]]$points


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
  print(of)
  pb$tick()
}

nu_opt_init=as.numeric(grid2[which.min(of_h_nu),])[1]
h_opt_init=as.numeric(grid2[which.min(of_h_nu),])[2]
pb$terminate()


c=seq(0.25,10,l=30)

grid=expand.grid(Nu,c)

of_c_nu=c()
pb <- progress_bar$new(total=nrow(grid))

for (i in 1:nrow(grid)){
  
  nu=as.numeric(grid[i,])[1]
  c=as.numeric(grid[i,])[2]
  of_c_nu=append(of_c_nu,
                 OF_c(y,X,c,nu,nu_init=nu_opt_init,h_init=h_opt_init)$logS)
  pb$tick()
}

nu_opt=as.numeric(grid[which.min(of_c_nu),])[1]
c_opt=as.numeric(grid[which.min(of_c_nu),])[2]
pb$terminate()

smooth_CDF=smooth_IDR_CDF(y,X,h=h_init,nu=nu_init,y_test=y_test,x_test=x_test)
smooth_CDF_local=smooth_IDR_CDF_h_opt(y,X,h_init=h_opt_init,nu_init=nu_opt_init,
                                      c=c_opt,nu=nu_opt,
                                      y_test=y_test,x_test=x_test)
smooth_CDF_local_uncorrected=smooth_IDR_CDF_h_opt(y,X,h_init=h_opt_init,nu_init=nu_opt_init,
                                                  c=c_opt,nu=nu_opt,
                                                  y_test=y_test,x_test=x_test,
                                                  increasing = FALSE)

plot(pred)
lines(y_test,smooth_CDF,col=2)
lines(y_test,smooth_CDF_local,col=3)
lines(y_test,smooth_CDF_local_uncorrected,col=3,lty=2)

legend("bottomleft",c("IDR", "smooth IDR global h", "smooth IDR local h",
                      "smooth IDR local h, uncorrected"),
       col=c("blue",2,3,3),lty=c(1,1,1,2))




