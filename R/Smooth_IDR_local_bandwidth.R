library(isodistrreg)
library(progress)

# In this file, we implement the SMOOTH IDR CDF AND DENSITY with automated local bandwidth selection as intruduced in the summary. 
# this requires selection of initial values of nu and h and c. These can be chosen a priori or by minimizing OF. 
# Introducing a local bandwidth has as consequence, that the resulting density does no longer integrate to one and that the resulting CDF
# is no longer guarenteed to be increasing. For the latter, we provide an adaptation that corrects this behaviour. For the former,
# it could be enough to normalize by dividing the estimates of the density with the sum of all estimated densities.
#

smooth_IDR_CDF_h_opt=function(y,x,y_test=y,x_test,c=1,
                              h_init=(log(length(x))/length(x))^(1/10),
                              nu=2.5,nu_init=nu,progress=FALSE,
                              increasing=TRUE){
 
  fit_idr=idr(y,data.frame(x))
  X=x
  unique_order_y=sort(unique(y))
  x=sort(unique(x))
  n=length(x)
  if (nu!=Inf){
    K_h=function(u,h){1/h*(1/(nu*pi)^(0.5))*gamma((nu+1)/2)/gamma(nu/2)*(1+(u/h)^2/nu)^(-(nu+1)/2)}
  }
  else{
    K_h=function(u,h){1/h*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)}
    }
  if (nu_init!=Inf){
    K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
    (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2)
  )}}
  else{
    K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-(u/h)^2/2)*(u^2/h^2-1)}
  }
  
  
  
  #compute w_j(x_test)
  w=NULL
  if (x_test>min(x) & x_test<max(x)){
    for(i in 1:(n-1)){
      if (x[i]<x_test& x_test<x[i+1]){
        
        w=c(1/(x[i+1]-x[i])*((x_test-x[i])*fit_idr$cdf[i+1,1]+
                               (x[i+1]-x_test)*fit_idr$cdf[i,1]))
        
        for (j in 2:length(unique_order_y)){
          w_j=1/(x[i+1]-x[i])*((x_test-x[i])*fit_idr$cdf[i+1,j]+
                                 (x[i+1]-x_test)*fit_idr$cdf[i,j])
          w=append(w,w_j)
          
        }
      }else if (x[i]==x_test){
        w=c(fit_idr$cdf[i,1])
        for (j in 2:length(unique_order_y)){
          w_j=fit_idr$cdf[i,j]
          w=append(w,w_j)
        }
        
      }
      
    }}else{if(x_test<=min(x)){
      i=1
      w=c(fit_idr$cdf[i,1])
      for (j in 2:length(unique_order_y)){
        w_j=fit_idr$cdf[i,j]
        w=append(w,w_j)
      }
    }
      else{ if (x_test>=max(x)){
        i=n
        w=c(fit_idr$cdf[i,1])
        for (j in 2:length(unique_order_y)){
          w_j=fit_idr$cdf[i,j]
          w=append(w,w_j)
        }
      }
      }}
  
  smooth_IDR_x_test=c()
  if(progress){pb <- progress_bar$new(total=length(y_test))}
  
  for (i in 1:length(y_test)){
    Y=sort(y_test)[i]
    K2=c()
    for (yy in unique_order_y){
      K2=append(K2,K_h_2(yy-Y,h_init))}
    
    
    second_der= w%*%K2
    kappa0=K_h(0,1)
    eps_n=(log(n)/n)^1/3
    if(nu!=Inf){var_gamma=nu/(nu-2)}else{var_gamma=1}
    
    h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var_gamma))^(1/3)
    K=c()
   # print(c(i,second_der,h_opt_Y))
    for(j in 1:(length(unique_order_y)-1)){
      
      #K=append(K,integrate(function(u){K_h(Y-u,h_opt_Y)},lower = unique_order_y[j],
       #                    upper = unique_order_y[j+1],stop.on.error = FALSE)$value)
      K=append(K,ifelse(nu!=Inf,integrate(function(u){K_h(Y-u,h_opt_Y)},
                                          lower = unique_order_y[j],
                                          upper = unique_order_y[j+1],
                                          stop.on.error = FALSE)$value,
                        0.5*(erf((Y-unique_order_y[j])/h_opt_Y)-
                                  erf((Y-unique_order_y[j+1])/h_opt_Y))))
      #K=append(K,mean(c(K_h(Y-unique_order_y[j+1]),K_h(Y-unique_order_y[j]))))
    }
    
    #K=append(K,integrate(function(u){K_h(Y-u,h_opt_Y)},lower = unique_order_y[length(unique_order_y)],
     #                    upper =Inf,stop.on.error = FALSE)$value)
    K=append(K,ifelse(nu!=Inf,integrate(function(u){K_h(Y-u,h_opt_Y)},
                                lower = unique_order_y[length(unique_order_y)],
                                upper = Inf)$value,
                      0.5*(erf((Y-unique_order_y[length(unique_order_y)])/h_opt_Y)+1)))
    #K=append(K,K_h(Y-unique_order_y[length(unique_order_y)]))
    #print(h_opt_Y)
    smooth_IDR_x_test=append(smooth_IDR_x_test,
                             w%*%K)
    
    if(progress){pb$tick()}
    }
  
  
  if(increasing){
    y_test=sort(y_test)
    u=smooth_IDR_x_test
    #u[length(u)]=ifelse(u[length(u)-1]>=0.995,1,u[length(u)-1])
    count=0
    
    
    while ((sum(u==sort(u))!=length(u))&count<=100){
      for (i in 2:length(u)){
        if (u[i]<u[i-1]){
         
          k=0
          for (j in 1:(i-1)){
            if(u[j]<u[j+1]){k=j}
          }
          if(k==0){
            l=min(which(u[2:length(u)]-u[1:(length(u)-1)]>0))+1
            m=u[l]/(y_test[l]-y_test[1])
            u[1]=0
            for (s in 2:(l-1)){
              u[s]=(y_test[s]-y_test[1])*m
            }
          }else{
            l=min(which(u[i-1]-u[(i-1):length(u)]>0))
            o=ifelse(i-1+l<=length(u),i-1+l,length(u))
            m=(u[o]-u[i-1])/(y_test[o]-y_test[i-1])
            for(s in i:(o-1)){
              u[s]=m*(y_test[s]-y_test[i-1])+u[i-1]
            }
          }
        }
      }
      count=count+1
      
    }
    smooth_IDR_x_test=u
    #print(count)
   
    }
    
  return(smooth_IDR_x_test)
  
} 



smooth_IDR_density_h_opt=function(y,x,x_test,y_test=y,
                                  c=1,h_init=0.5,nu=2.5,
                                  nu_init=nu,normalize=FALSE){
  fit_idr=idr(y,data.frame(x))
  n=length(unique(x))
  x=sort(unique(x))
  
  unique_order_y=sort(unique(y))
  
  if(nu!=Inf){
    K_h=function(u,h){1/h*(1/(nu*pi)^(0.5))*gamma((nu+1)/2)/gamma(nu/2)*(1+(u/h)^2/nu)^(-(nu+1)/2)}
    }else{
    K_h=function(u,h){1/h*(2*pi)^(-0.5)*exp(-(u/h)^2/2)}
   }
  
  if(nu_init!=Inf){
     K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
      (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
  }else{
     K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-(u/h)^2/2)*(u^2/h^2-1)}
  }
  
  
  #if(nu_init!=Inf){
  #  K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
  #    (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
  #}else{
  #  K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)*(u^2/h^2-1)}
  #}
  #compute w_j(x_test)
  w=NULL
  if (x_test>min(x) & x_test<max(x)){
    for(i in 1:(n-1)){
      if (x[i]<x_test& x_test<x[i+1]){
        w=c(1/(x[i+1]-x[i])*((x_test-x[i])*fit_idr$cdf[i+1,1]+
                               (x[i+1]-x_test)*fit_idr$cdf[i,1]))
        for (j in 2:length(unique_order_y)){
          w_j=1/(x[i+1]-x[i])*((x_test-x[i])*fit_idr$cdf[i+1,j]+
                                 (x[i+1]-x_test)*fit_idr$cdf[i,j])-
            1/(x[i+1]-x[i])*((x_test-x[i])*fit_idr$cdf[i+1,j-1]+
                               (x[i+1]-x_test)*fit_idr$cdf[i,j-1])
          w=append(w,w_j)
          
        }
      }else if (x[i]==x_test){
        w=c(fit_idr$cdf[i,1])
        for (j in 2:length(unique_order_y)){
          w_j=fit_idr$cdf[i,j]-fit_idr$cdf[i,j-1]
          w=append(w,w_j)
        }
        
      }
      
    }}else{if(x_test<=min(x)){
      i=1
      w=c(fit_idr$cdf[i,1])
      for (j in 2:length(unique_order_y)){
        w_j=fit_idr$cdf[i,j]-fit_idr$cdf[i,j-1]
        w=append(w,w_j)
      }
    }
      else{ if (x_test>=max(x)){
        i=n
        w=c(fit_idr$cdf[i,1])
        for (j in 2:length(unique_order_y)){
          w_j=fit_idr$cdf[i,j]-fit_idr$cdf[i,j-1]
          w=append(w,w_j)
        }
      }
      }}
  
  smooth_IDR_x_test=c()
  for (Y in sort(y_test)){
    #second_der=second_derivative_f(y,x,h_init,Y,x_test)
    K2=sapply(sort(unique(y)),function(u){K_h_2(u-Y,h_init)})
    second_der= w%*%K2
    kappa0=K_h(0,1)
    eps_n=(log(n)/n)^1/3
    if(nu!=Inf){var_gamma=nu/(nu-2)}else{var_gamma=1}
    
    h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var_gamma))^(1/3)
    
    smooth_IDR_x_test=append(smooth_IDR_x_test,
                             w%*%sapply(unique_order_y,
                                        function(u){K_h(Y-u,h_opt_Y)}
                                        
                             ))}
  if(normalize){
    integrand=function(s){
      K2=sapply(sort(unique(y)),function(u){K_h_2(u-s,h_init)})
      
      second_der= w%*%K2
      kappa0=K_h(0,1)
      eps_n=(log(n)/n)^1/3
      if(nu!=Inf){var_gamma=nu/(nu-2)}else{var_gamma=1}
      
      h_opt_Y=c*(4*kappa0*eps_n/(abs(second_der)*var_gamma))^(1/3)
      return(w%*%sapply(unique_order_y,function(u){K_h(s-u,h_opt_Y)}))
    }
    range=range(y)+c(-100,100)
    int=integrate(Vectorize(integrand),lower=range[1],upper=range[2])$value
    
  }
  return(list(density=ifelse(normalize,1/int,1)*smooth_IDR_x_test
              ,weights=w,integral=ifelse(normalize,int,NA)))
  
}
