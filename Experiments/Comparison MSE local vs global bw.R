#install.packages("devtools")
library(devtools)
#loads the neccessary functions
source_url("https://raw.githubusercontent.com/timst91/Smooth-IDR/main/R/all.R")

id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


H_init=seq(0.01,2,l=5)
C=c(0.25,0.5,1,2,4)
NU=c(2.01,3,5,10,Inf)
H=seq(0.01,4,l=10)


grid_local=expand.grid(NU,C,H_init)
grid_global=expand.grid(NU,H)

y_test=seq(0,80,l=50)
x_test=6

N=seq(50,1000,by=150)

CDFs_local=matrix(ncol=length(N),nrow=length(y_test))

CDFs_global_CV=matrix(ncol=length(N),nrow=length(y_test))

CDFs_global_OF=matrix(ncol=length(N),nrow=length(y_test))


for(i in 1:length(N)){
  n=N[i]
  X=runif(n,0,10)
  y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))
  
  
  CV_local=apply(grid_local,1,function(w){
    CV_c(y=y,x=X,c=w[2],nu=w[1],h_init=w[3],progress = 0)
  })
 
  CV_global=apply(grid_global,1,function(w){
    CV_h(y=y,x=X,nu=w[1],h=w[2],progress = 0)
  })
  
  OF=apply(grid_global,1,function(w){
    OF_h(y=y,x=X,nu=w[1],h=w[2],progress = 0)
  })

  par_local=as.numeric(grid_local[which.min(CV_local),])
  par_global_cv=as.numeric(grid_global[which.min(CV_global),])
  par_global_of=as.numeric(grid_global[which.min(OF),])
  
  
  CDFs_local[,i]=smooth_IDR_CDF_h_opt(y,X,x_test = x_test,y_test = y_test,
                                      c=par_local[2],nu=par_local[1],
                                      h_init = par_local[3])$cdf
  
  CDFs_global_OF[,i]=smooth_IDR_CDF(y,X,x_test = x_test,y_test = y_test,
                                    h=par_global_of[2],nu=par_global_of[1],
  )
  CDFs_global_CV[,i]=smooth_IDR_CDF(y,X,x_test = x_test,y_test = y_test,
                                    h=par_global_cv[2],nu=par_global_cv[1],
  )
  
}


CDFs=list(CDFs_local=CDFs_local,CDFs_global_OF=CDFs_global_OF,
          CDFs_global_CV=CDFs_global_CV)

save(
  list = "CDFs",
  # save in new file for each id
  file = paste0("CDFs_MSE", id, ".rda")
)




#############Evaluation######################



res <- vector("list", 100) #100 simulations on Euler
for (id in 1:100) {
  load(paste0("CDFs_MSE", id, ".rda"))
  res[[id]] <- CDFs
}
res <- do.call(rbind, res)

ids=1:100

CDFs_local=0
for (i in 1:length(ids)){
  CDFs_local=CDFs_local+1/length(ids)*res[[i]]
}

CDFs_of=0
for (i in (length(ids)+1):(2*length(ids))){
  CDFs_of=CDFs_of+1/length(ids)*res[[i]]
}

CDFs_cv=0
for (i in (2*length(ids)+1):(3*length(ids))){
  CDFs_cv=CDFs_cv+1/length(ids)*res[[i]]
}



y_test=seq(0,80,l=50)
N=seq(50,1000,by=150)

#ind=1:50  

ind_N=c(1,2,4,5,7)  #less values for N for visibility
cols=c(2,3,4,5,6)
#ind_N=1:length(N)
#cols=2:8
labels=sapply(ind_N, function(i) substitute(paste("n", " = ", n), list( n= N[i])))



par(mfrow=c(1,3))
plot(y_test[ind],pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8))[ind],type="l",
     ylab="",xlab="Threshold")

legend("right",legend=labels,
       col=cols,
       lty=1,bty="n",ncol = 1,seg.len = 1,text.width = 5)
#ncol=2

for(i in 1:length(ind_N)){
  j=ind_N[i]
  lines(y_test[ind],CDFs_local[ind,j],col=cols[i],lwd=.5)
}

plot(y_test[ind],pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8))[ind],type="l",
     ylab="",xlab="Threshold")


for(i in 1:length(ind_N)){
  j=ind_N[i]
  lines(y_test[ind],CDFs_of[ind,j],col=cols[i])
}

plot(y_test[ind],pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8))[ind],type="l",
     ylab="",xlab="Threshold")


for(i in 1:length(ind_N)){
  j=ind_N[i]
  lines(y_test[ind],CDFs_cv[ind,j],col=cols[i])
}

par(mfrow=c(1,1))

########## MSE plot ###############

MSE_local_N=apply(CDFs_local[,ind_N],2,function(x){
  sum((x-pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8)))^2)})

MSE_of_N=apply(CDFs_of[,ind_N],2,function(x){
  sum((x-pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8)))^2)})


MSE_cv_N=apply(CDFs_cv[,ind_N],2,function(x){
  sum((x-pgamma(y_test,shape=sqrt(x_test),scale=min(max(x_test,2),8)))^2)})

plot(N[ind_N],MSE_local_N,type="b",pch=3,ylim=c(0.3,0.6))
lines(N[ind_N],MSE_cv_N,type="b",pch=2,col=2)

lines(N[ind_N],MSE_of_N,type="b",pch=4,col=3)

