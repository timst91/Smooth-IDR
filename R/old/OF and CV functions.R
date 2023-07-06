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


# Function that computes OF in terms of c for smooth IDR with optimal local bandwidth selection procedure

OF_c=function(y,x,c,nu,nu_init=2.5,h_init=(log(length(unique(x)))/length(unique(x)))^(1/10),progress=FALSE){
  sort_ind= match(sort(unique(x)), x)
  n=length(sort_ind)
  x=x[sort_ind]
  y=y[sort_ind]
  idr_fit=idr(y,data.frame(x))
  if(progress){pb <- progress_bar$new(total=n)}
  
  if(nu!=Inf){
    K_h=function(u,h){1/h*(1/(nu*pi)^(0.5))*gamma((nu+1)/2)/gamma(nu/2)*(1+(u/h)^2/nu)^(-(nu+1)/2)}
  }else{
    K_h=function(u,h){1/h*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)}
    }
  if(nu_init!=Inf){
    K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
    (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
  }else{
      K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)*(u^2/h^2-1)}
      }

 
  cdf=idr_fit$cdf
  count=0
  logS=0
  for (i in 1:n){
    w=c(cdf[i,1])
    for (j in 2:length(unique(y))){
      w_j=cdf[i,j]-cdf[i,j-1]
      w=append(w,w_j)
    }
    
    K2=sapply(sort(unique(y)),function(u){K_h_2(u-y[i],h_init)})
    
    second_der= w%*%K2
    kappa0=K_h(0,1)
    
    eps_n=(log(n)/n)^1/3
    if(nu!=Inf){var_gamma=nu/(nu-2)}else{var_gamma=1}
    
    h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var_gamma))^(1/3)
    
    K=sapply(sort(unique(y)),function(u){K_h(u-y[i],h_opt_Y)})
    
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
    
   if(sum(w_i%*%K) !=0){
      logS=logS-1/n*
        log(w_i%*%K)
    }else{count=count+1}
    
    if(progress){pb$tick()}
  }
  #print(logS)
  return(n/(n-count)*logS)
  
  
}

CV_c=function(y,x,c){
  sort_ind= match(sort(unique(x)), x)
  n=length(sort_ind)
  x_unique=x[sort_ind]
  logS=0
  pb <- progress_bar$new(total=length(x_unique))
  
  
  for (i in 1:length(x_unique)){
    model=smooth_IDR_density_h_opt(y=y[sort_ind][-i],x=x_unique[-i],c=c,
                               y_test = y[sort_ind][i],x_test =x_unique[i])
    logS=logS-1/n*log(model)
    pb$tick()
  }
  # Close the progress bar
  
  return(logS)
}
