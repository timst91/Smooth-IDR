y_test=seq(0,60,l=100)

N=seq(100,4000,l=20)
MSE=c()
MSE_opt=c()
x_test=5
real_cdf=pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8))
nsim=5

H=seq(0.1,2,l=10)
Nu=c(seq(2.01,10,l=10),Inf)

C=c(0.25,0.5,1,2,4)
grid=expand.grid(Nu,C)
grid2=expand.grid(Nu,H)

for (n in N2){
  mse=0
  mse_opt=0
  for (j in 1:nsim){
    y=c()
    X=sort(runif(n,0,10))
    
    for (x in X) {y=append(y,rgamma(1,shape=sqrt(x),scale=min(max(x,2),8)))}
    
    of_h_nu=c()
    
    for (i in 1:nrow(grid2)){
      
      nu=as.numeric(grid2[i,])[1]
      h=as.numeric(grid2[i,])[2]
      of=OF_h(y,X,h,nu)
      of_h_nu=append(of_h_nu,
                     of)
    }
    nu_opt_init=as.numeric(grid2[which.min(of_h_nu),])[1]
    h_opt_init=as.numeric(grid2[which.min(of_h_nu),])[2]
    
    print(c("nu_init=",nu_opt_init,"h_init=",h_opt_init))
    smooth.idr.cdf=smooth_IDR_CDF(y,X,h=h_opt_init,nu=nu_opt_init, y_test = y_test,x_test=x_test)
    
    of_c_nu=c()
    
    for(i in 1:nrow(grid)){
      nu=as.numeric(grid[i,])[1]
      c=as.numeric(grid[i,])[2]
      of=OF_c(y,X,c,nu,nu_init = nu_opt_init,h_init = h_opt_init)
      of_c_nu=append(of_c_nu,of)
      
    }
    nu_opt=as.numeric(grid[which.min(of_c_nu),])[1]
    c_opt=as.numeric(grid[which.min(of_c_nu),])[2]
    print(c("nu=",nu_opt,"c=",c_opt))
    
    
    smooth.idr.cdf.h.opt=smooth_IDR_CDF_h_opt(y=y,x=X,c=c_opt, h_init = h_opt_init,
                                              nu=nu_opt,nu_init = nu_opt_init,
                                              y_test = y_test,x_test=x_test)
    mse=mse+1/nsim*mean((smooth.idr.cdf-real_cdf)^2)
    mse_opt=mse_opt+1/nsim*mean((real_cdf-smooth.idr.cdf.h.opt)^2)
  }
  #lines(y_test,smooth.idr.cdf)
  print(mse)
  MSE=append(MSE,mse)
  MSE_opt=append(MSE_opt,mse_opt)
  
}

plot(N,MSE,type='l')
lines(N,MSE_opt,col=2)
legend("topright",legend=c("local bw","global bw"),col=c(2,1),lty=1)

