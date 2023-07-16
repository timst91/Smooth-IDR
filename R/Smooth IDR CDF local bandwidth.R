library(cubature)
library(isodistrreg)

smooth_IDR_CDF_h_opt=function(y,x,x_test,
                              y_test=y,
                              c=1,
                              h_init=(log(length(x))/length(x))^(1/9),
                              nu=2.5,
                              nu_init=nu,
                              rel_tol=1e-3,
                              progress=FALSE,tol=1e-4,
                              lower="automatic"){
  
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
  y_test=sort(unique(y_test))
  m_test=length(y_test)
  
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
  
  if(progress){pb=progress_bar$new(total=m_test+1)}
  
  int=smooth_IDR_density_h_opt(y,x,x_test,y_test=mean(y),
                               c=c,
                               nu=nu,
                               nu_init = nu_init,
                               h_init=h_init,
                               normalize = TRUE)$integral
  
  if(progress){pb$tick()}
  
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
  
  
  kappa0=K_h(0,1)
  eps_n=(log(n)/n)^1/3
  
  if(nu!=Inf){var=nu/(nu-2)}else{var=1}
  
  
  h_opt=function(y){
    second_der= w%*%sapply(unique_order_y,function(u){K_h_2(u-y,h_init)})
    
    return(c*(4*kappa0*eps_n/(abs(second_der)*var))^(1/3))
  }
  
  
  if(lower=="automatic"){
    lower=0.5*min(c(y_test,y))-0.5*(max(c(y_test,y)))
  }
  
  
  
  
  nonzero=which(w!=0)
  w0=w[nonzero]
  unique_order_y0=unique_order_y[nonzero]
  cdf=numeric(m_test)+sum(w0*sapply(unique_order_y0,
                                  function(u){
                                    adaptIntegrate(
                                      Vectorize(function(w){K_h(w-u,h_opt(w))}),
                                      lowerLimit=lower,
                                      upperLimit=y_test[1],tol=tol,
                                      absError = tol)$integral
                                  }))
  
  if(progress){pb$tick()}
  
  if(m_test>=2){
    if(m_test==2){
      cdf[2]=cdf[1]+sum(w0*sapply(unique_order_y0,
                                  function(u){
                                    adaptIntegrate(
                                      Vectorize(function(w){K_h(w-u,h_opt(w))}),
                                      lowerLimit=y_test[1],
                                      upperLimit=y_test[2],tol=tol,
                                      absError = tol)$integral
                                  }))
      if(progress){pb$tick()}
    }else{
      for(i in 2:(m_test)){
       
        cdf[i]=cdf[i-1]+sum(w0*sapply(unique_order_y0,
                                      function(u){
                                        adaptIntegrate(
                                          Vectorize(function(w){K_h(w-u,h_opt(w))}),
                                          lowerLimit=y_test[i-1],
                                          upperLimit=y_test[i],tol=tol,
                                          absError = tol)$integral
                                      }))
        
        if(progress){pb$tick()}
      }
    }
  }
  
  if(progress){pb$terminate()}
  
  cdf=1/int*cdf
  return(list(cdf=cdf,integral=int))
  
}
