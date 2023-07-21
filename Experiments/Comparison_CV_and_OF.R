
#install.packages("devtools")
library(devtools)
#loads the neccessary functions
source_url("https://raw.githubusercontent.com/timst91/Smooth-IDR/main/R/all.R")




id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

N=seq(100,2000,l=20)
ln=length(N)

H=numeric(6)
H[2:6]=c(0.25,0.5,1,2,4)

of_h=matrix(nrow=ln,ncol=length(H))
cv_h=matrix(nrow=ln,ncol=length(H))
nu=Inf

time_of=numeric(ln)
time_cv=numeric(ln)

pb=progress_bar$new(total=ln*6)

for(i in 1:ln){
  
  n=N[i]
  H[1]=(log(n)/n)^(1/9)
  X=runif(n,0,10)
  y=rgamma(n,shape=sqrt(X),scale=min(max(X,2),8))
  
  for(j in 1:length(H)){
    
    h=H[j]
    
    of=OF_h(y,X,h,nu,time=TRUE)
    of_h[i,j]=of$logS
    time_of[i]=time_of[i]+1/length(H)*of$time
    
    cv=CV_h(y,X,h,nu,time=TRUE)
    cv_h[i,j]=cv$logS
    time_cv[i]=time_cv[i]+1/length(H)*cv$time
    pb$tick()
    if(j==1||j==6){
      cat(paste('j:',j,"| ", "of_logS:",of$logS,"| ",
                
              "cv_logS:",cv$logS,"| ",  ' time elapsed:',of$time,cv$time))
    }
    
  }
  
}
pb$terminate()

results=list(cv=cv_h,of=of_h,time=time)

save(
  list = "results",
  # save in new file for each id
  file = paste0("of_cv", id, ".rda")
)





############## Evaluation: #######################


res <- vector("list", 100)   #100 simulations on Euler

for (id in (1:100)) {
  load(paste0("of_cv", id, ".rda"))
  res[[id]] <- results
}
res <- do.call(rbind, res)

save(list = "res", file = "of_cv_results.rda")



mean_of_matrix=matrix(0,ncol=5,nrow=20)
mean_cv_matrix=mean_of_matrix


for(i in 1:100){
  
  mean_of_matrix=mean_of_matrix+1/nrow(mean_of_matrix)*res$of[[i]]
  mean_cv_matrix=mean_cv_matrix+1/nrow(mean_cv_matrix)*res$cv[[i]]
}



plot(N,apply(log(abs(mean_of_matrix-mean_cv_matrix)[,-1]),1,mean),
     ylim=c(-0.5,2.5),type='l',xlab="n",ylab="log(absolute difference)",col=4)

lines(N,mean(log(abs(mean_of_matrix-mean_cv_matrix)[,1])))



###################### SINGLE PLOT FOR VISUALIZATION OF BANDWIDTH CHOICE ################################
Hh=c(0.25,0.5,1:15)

ind_n=c(1,5,10)
ind_n=1
nsim=100

par(mfrow=c(1,length(ind_n)))
for( i in 1:length(ind_n)){
  n=N[ind_n[i]]
  of=0
  cv=0
  for(sim in 1:nsim){
    print(sim)
    X=runif(n,0,10)
    y=rgamma(n,shape=sqrt(X),scale = min(max(X,2),8))
    of=of+1/nsim*sapply(Hh,function(u){OF_h(y,X,u,Inf,progress = 1)})
    cv=cv+1/nsim*sapply(Hh,function(u){CV_h(y,X,u,Inf,progress = 1)})
  }
  
  plot(Hh,log(of),type="b",xlim=c(0,10.4),ylab="log(score)",xlab="h",col=2)
  legend("topright",c("OF","CV"),col=c(2,1),lty=1,bty="n")
  lines(Hh,log(cv),type="b",col=1) 
  abline(v=Hh[which.min(of)],lty=2,col=2)
  abline(v=Hh[which.min(cv)],col=1,lty=3)
}






