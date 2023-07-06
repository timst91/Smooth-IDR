smooth_IDR_CDF_h_opt=function(y,x,x_test,
                              y_test=y,
                              c=1,
                              h_init=(log(length(x))/length(x))^(1/10),
                              nu=2.5,
                              nu_init=nu,
                              max_h=1e2,
                              progress=FALSE
                              ){
  
  sort_ind=match(sort(unique(x)),x)
  
  x=sort(unique(x))
  y=y[sort_ind]
  
  unique_order_y=sort(unique(y))
  n=length(x)
  m=length(unique_order_y)
  
  fit.idr=idr(y,data.frame(x))
  pred.idr=predict(fit.idr,data=data.frame(x=x_test))
  
  # 2nd derivatives of kernel functions
  
  if (nu_init!=Inf){
    K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*
        gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*
        ((1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*
            (1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
    }else{
    K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-(u/h)^2/2)*(u^2/h^2-1)}
  }
  
  
  if (nu==Inf){
    
    K_h=function(x,h){dnorm(x,0,h)}
  }else{
    K_h=function(x,h){dt(x/h,nu)/h}
  }
  
  
  #compute cdfs
  
  
  fit.idr=idr(y,data.frame(x))
  pred.idr=predict(fit.idr,data=data.frame(x=x_test),digits=6)
  
  
  points=pred.idr[[1]]$points
  cdf.idr=pred.idr[[1]]$cdf
  
  cdf=numeric(m)
  
  for(j in 1:m){
    y=unique_order_y[j]
    ind_y=ifelse(points[1]<=y,max(which(points <= y)),0)
    
    cdf[j]=ifelse(ind_y>0,cdf.idr[ind_y],0)
    
  }
  
  #compute w_j(x_test)
  
 
  w=numeric(m)
  
  for(j in 1:m){
    y=unique_order_y[j]
    ind_y=ifelse(points[1]<=y,max(which(points <= y)),0)
    
    if(j==1){
      w[j]=ifelse(ind_y>0,cdf.idr[ind_y],0)
    }else{
      
      y1= unique_order_y[j-1]
      
      ind_y_1=ifelse(points[1]<=y1,max(which(points <= y1)),0)
      
      
      if(ind_y>0){
        if(ind_y_1>0){
          w[j]=ifelse(ind_y!=ind_y_1,cdf.idr[ind_y]-cdf.idr[ind_y_1],0)
        }else{
          w[j]=cdf.idr[ind_y]}
      }else{w[j]=0}
      
    }
    
    
  }
  
  K_int=numeric(m)
  
  
  y_test=sort(y_test)
  m_test=length(y_test)
  
  smooth_IDR_x_test=numeric(m_test)
  if(progress){pb = progress_bar$new(total=m_test)}
  u=numeric(m_test)
  
  kappa0=K_h(0,1)
  eps_n=(log(n)/n)^1/3
  
  if(nu!=Inf){var=nu/(nu-2)}else{var=1}
    
  for (i in 1:m_test){
    Y=y_test[i]
   
    second_der= as.numeric(w%*%sapply(unique_order_y,function(u){K_h_2(Y-u,h_init)}))
    
   
    h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var))^(1/3)
    h_opt_Y=ifelse(max_h>h_opt_Y,h_opt_Y,max_h)
    u[i]=h_opt_Y
    
    if(nu== Inf){
      
      for(j in 1:(m-1)){
        K_int[j]=pnorm(unique_order_y[j+1],Y,h_opt_Y)-
          pnorm(unique_order_y[j],Y,h_opt_Y)
      }
      K_int[m]=1-pnorm(unique_order_y[m],Y,h_opt_Y)
    }else{
      
      for(j in 1:(m-1)){
        K_int[j]=pt((unique_order_y[j+1]-Y)/h_opt_Y,nu)-
          pt((unique_order_y[j]-Y)/h_opt_Y,nu)
      }
      
      K_int[m]=1-pt((unique_order_y[m]-Y)/h_opt_Y,nu)
    }
    if(progress){pb$tick()}
    
    smooth_IDR_x_test[i]=cdf%*%K_int
    
  }

  if(progress){pb$terminate()}
  return(list(cdf=smooth_IDR_x_test,opt_bw=u))
  
} 


