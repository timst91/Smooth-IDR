library(tictoc)
library(isodistrreg)

OF_h=function(y,x,h,nu,time=FALSE,progress=FALSE){
  if(time){tic()}
  
  sort_ind=match(sort(unique(x)),x)
  
  x=sort(unique(x))
  y=y[sort_ind]
  
  unique_order_y=sort(unique(y))
  n=length(x)
  m=length(unique_order_y)
  
  fit.idr=idr(y,data.frame(x))
  cdf=fit.idr$cdf
  
  if(progress){pb = progress_bar$new(total=n)}
  
  if (nu==Inf){
    K_h=function(x){dnorm(x,0,h)}
  }else{
    K_h=function(x){dt(x/h,nu)/h}
  }
  
  
 
  logS=numeric(n)
  count=0
  wk=c()
  
  for (i in 1:n){
    
    w_1=cdf[i,c(1,1:(m-1))]
    w_1[1]=0
    w_i=cdf[i,]-w_1
    j_i=which(match(unique_order_y, y[i])==1)
    w_j_i=w_i[j_i]
    w_i[j_i]=0
    if(w_j_i!=1){
      w_i=w_i/(1-w_j_i)
    }else{
      count=count+1
      w_i=w_i}
    
    
    K=sapply(unique_order_y,function(u){K_h(u-y[i])})
    wik=w_i%*%K
    
    if(sum(wik)!=0){
      logS[i]=-log(wik)
      wk=append(wk,wik)
      
    }else{
      count=count+1
      if(length(wk)>=1){
        wik=ifelse(length(wk)<=5,mean(wk),
                   mean(wk[(length(wk)-5):length(wk)]))
        
        logS[i]=-log(wik)
        wk=append(wk,wik)
      }
      
      
    }
    
    if(progress){pb$tick()}
    
  }
  
  if(progress){pb$terminate()}
  if(time){Time=toc(quiet = TRUE)
  
  
  return(list(logS=n/(n-count)*mean(logS),time=unname(Time$toc-Time$tic)))
  }else{
    return(mean(logS))
  }
}


CV_h=function(y,x,h,nu,time=FALSE,progress=FALSE){
  
  if(time){tic()}
  
  sort_ind=match(sort(unique(x)),x)
  
  x=sort(unique(x))
  y=y[sort_ind]
  
  unique_order_y=sort(unique(y))
  n=length(x)
  
  logS=numeric(n)
  
  if (progress){pb = progress_bar$new(total=n)}
  
  
  for (i in 1:n){
    model=smooth_IDR_density(y=y[-i],x=x[-i],
                             h=h,nu=nu,
                             y_test = y[i],
                             x_test =x[i])
    logS[i]=-log(model)
    
    if(progress){pb$tick()}
  }
  if(progress){pb$terminate()}
  if(time){
    time=toc(quiet=TRUE)
    return(list(logS=mean(logS),time=unname(time$toc-time$tic)))
  }else{
    return(mean(logS))
  }
  
  
}


# Function that computes OF in terms of c for smooth IDR with optimal local bandwidth selection procedure

OF_c=function(y,x,c,nu,
              nu_init=2.5,
              h_init=(log(length(unique(x)))/length(unique(x)))^(1/10),
              time=FALSE,progress=FALSE){
  
  if(time){tic()}
  
  sort_ind=match(sort(unique(x)),x)
  
  x=sort(unique(x))
  y=y[sort_ind]
  
  unique_order_y=sort(unique(y))
  n=length(x)
  m=length(unique_order_y)
  
  fit.idr=idr(y,data.frame(x))
  cdf=fit.idr$cdf
  
  if(progress){pb <- progress_bar$new(total=n)}
  
  if (nu==Inf){
    K_h=function(x){dnorm(x,0,h)}
  }else{
    K_h=function(x){dt(x/h,nu)/h}
  }
  
  
  if(nu_init!=Inf){
    K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
      (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
  }else{
    K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)*(u^2/h^2-1)}
  }
  
  
  
  count=0
  logS=numeric(n)
  
  kappa0=K_h(0,1)
  
  eps_n=(log(n)/n)^1/3
  var=ifelse(nu!=Inf,nu/(nu-2),1)
  
  for (i in 1:n){
    w=cdf[i,1]+numeric(m)
    for (j in 2:m){
      w_j=cdf[i,j]-cdf[i,j-1]
      w=append(w,w_j)
    }
    
    K2=sapply(unique_order_y,function(u){K_h_2(u-y[i],h_init)})
    
    second_der= w%*%K2
    
    h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var))^(1/3)
    
    K=sapply(unique_order_y,function(u){K_h(u-y[i],h_opt_Y)})
    
    w_1=cdf[i,c(1,1:(length(unique(y))-1))]
    w_1[1]=0
    w_i=cdf[i,]-w_1
    j_i=which(match(unique_order_y, y[i])==1)
    w_j_i=w_i[j_i]
    w_i[j_i]=0
    if(w_j_i!=1){
      w_i=w_i/(1-w_j_i)
    }else{
      count=count+1
      w_i=w_i}
    
    if(sum(w_i%*%K) !=0){
      logS[i]=-log(w_i%*%K)
    }else{
      count=count+1
      logS[i]= logS[i-1]
      }
    
    if(progress){pb$tick()}
  }
  
  if(progress){pb$terminate()}
  
  if(time){
    time=toc(quiet=TRUE)
    return(list(logS=n/(n-count)*mean(logS),time=unname(time$toc-time$tic)))
  }else{
    return(mean(logS))
  }
  
}

CV_c=function(y,x,c,progress=FALSE){
  sort_ind= match(sort(unique(x)), x)
  x_unique=x[sort_ind]
  
  y=y[sort_ind]
  n=length(x_unique)
  
  logS=numeric(n)
  if(progress){pb= progress_bar$new(total=length(x_unique))}
  
  
  for (i in 1:n){
    model=smooth_IDR_density_h_opt(y=y[-i],x=x_unique[-i],c=c,
                                   y_test = y[sort_ind][i],
                                   x_test =x_unique[i])$density
    logS[i]=-log(model)
    if(progress){pb$tick()}
  }
  if(progress){pb$terminate()}
  
  return(mean(logS))
}
