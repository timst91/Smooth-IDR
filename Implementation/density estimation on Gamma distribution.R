y=c()
n=200
X=sort(runif(n,0,10))

for (x in X) {y=append(y,rgamma(1,shape=sqrt(x),scale=min(max(x,2),8)))}

x_test=5



idr_fit=idr(y,data.frame(X))

pred=predict(idr_fit,data=data.frame(X=x_test))

#y_test=pred[[1]]$points
y_test=seq(-5,37,l=100)

true_density=dgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8))


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

dens_IDR=smooth_IDR_density(y,X,h=h_opt_init,nu=nu_opt_init,y_test=y_test,x_test=x_test)
dens_IDR_local=smooth_IDR_density_h_opt(y,X,h_init=h_opt_init,
                                      c=c_opt,nu=nu_opt,
                                      y_test=y_test,x_test=x_test)


plot(y_test,dens_IDR,col=2,type='l',ylim=c(0,0.1))
lines(y_test,dens_IDR_local$density,col=1)
lines(y_test,true_density,lty=2,col="darkgreen")

Legend1= substitute(paste("IDR with global bw, ", h[opt], " = ", H, " , ", nu[opt]," = ",nU)
                    , list(H = h_opt_init, nU=nu_opt_init))

Legend2=substitute(paste("IDR with local bw, ", c[opt], " = ", C,
                         " , ", nu[opt]," = ",nU), list(C = c_opt, nU=nu_opt))



legend("topright",c("IDR", Legend1, Legend2,
                 expression(f[5](y))),
       col=c(2,1,"darkgreen"),lty=c(1,1,2))



