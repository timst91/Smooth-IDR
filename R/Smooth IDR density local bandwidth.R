library(isodistrreg)
library(cubature)

smooth_IDR_density_h_opt=function(y,x,x_test,
                                  y_test=y,
                                  c=1,
                                  h_init=0.5,
                                  nu=2.5,
                                  nu_init=nu,
                                  normalize=FALSE,
                                  progress=FALSE){
  
  sort_ind=match(sort(unique(x)),x)
  
  x=sort(unique(x))
  y=y[sort_ind]
  
  unique_order_y=sort(unique(y))
  n=length(x)
  m=length(unique_order_y)
  
  fit.idr=idr(y,data.frame(x))
  pred.idr=predict(fit.idr,data=data.frame(x=x_test))
  
  cdf=pred.idr[[1]]$cdf
  points=pred.idr[[1]]$points
  
  if (nu==Inf){
    K_h=function(x,h){dnorm(x,0,h)}
  }else{
    K_h=function(x,h){dt(x/h,nu)/h}
  }
  
  if(nu_init!=Inf){
    K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
      (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
  }else{
    K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-(u/h)^2/2)*(u^2/h^2-1)}
  }
  
  
  w=numeric(m)
  
  for(j in 1:m){
    
    y=unique_order_y[j]
    
    ind_y=ifelse(points[1]<=y,max(which(points <= y)),0)
    
    if(j==1){
      w[j]=ifelse(ind_y>0,cdf[ind_y],0)
    }else{
      
      y1= unique_order_y[j-1]
      
      ind_y_1=ifelse(points[1]<=y1,max(which(points <= y1)),0)
      
      if(ind_y>0){
        if(ind_y_1>0){
          w[j]=ifelse(ind_y!=ind_y_1,cdf[ind_y]-cdf[ind_y_1],0)
        }else{
          w[j]=cdf[ind_y]}
      }else{w[j]=0}
      
    }
    
    
  }
  nonzero=which(w!=0)
  w=w[nonzero]
  unique_order_y=unique_order_y[nonzero]
  
  K_int=numeric(m)
  
  y_test=sort(y_test)
  m_test=length(y_test)
  
  smooth_IDR_x_test=numeric(m_test)
  
  if(progress){pb = progress_bar$new(total=m_test)}
  
  kappa0=K_h(0,1)
  eps_n=(log(n)/n)^1/3
  
  if(nu!=Inf){var=nu/(nu-2)}else{var=1}
  
  for (i in 1:m_test){
    Y=y_test[i]
    
    K2=sapply(unique_order_y,function(u){K_h_2(u-Y,h_init)})
    second_der= w%*%K2
    
    h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var))^(1/3)
    
    smooth_IDR_x_test[i]=w%*%sapply(unique_order_y,
                                    function(u){K_h(Y-u,h_opt_Y)})
    
    if(progress){pb$tick()}
    
  }
  
  if(normalize){
    integrand=function(s){
      K2=sapply(unique_order_y,function(u){K_h_2(u-s,h_init)})
      
      second_der= w%*%K2
      h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var))^(1/3)
      
      return(w%*%sapply(unique_order_y,function(u){K_h(s-u,h_opt_Y)}))
    }
    
    
    range=range(unique_order_y)+c(-10,10)
    
    
    #int=integrate(Vectorize(integrand),lower=range[1],upper=range[2],
    #              subdivisions = ifelse(range_int=="finite",500,2000)
    #)$value
    int=adaptIntegrate(Vectorize(integrand),lowerLimit =range[1],upperLimit =range[2])$integral
    
  }
  
  if(progress){pb$terminate()}
  
  return(list(density=ifelse(normalize,1/int,1)*smooth_IDR_x_test,
              weights=w,
              integral=ifelse(normalize,int,NA)))
  
}
