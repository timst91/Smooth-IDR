library(progress)
library(pracma)
library(isodistrreg)
library(tictoc)


smooth_IDR_CDF=function(y,x,h=3*(log(length(x))/length(x))^(1/9),nu=2.5,y_test=y,x_test,progress=FALSE){
  #tic()
  fit_idr=idr(y,data.frame(x))
  
  unique_order_y=sort(unique(y))
  x=sort(unique(x))
  
  n=length(x)
  
  if (nu!=Inf){
    K_h=function(u){1/h*(1/(nu*pi)^(0.5))*gamma((nu+1)/2)/gamma(nu/2)*
        (1+(u/h)^2/nu)^(-(nu+1)/2)}
  }
  else{
    K_h=function(u){1/h*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)}
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
  if(progress){pb <- progress_bar$new(total=length(unique(y_test)))}
  
  smooth_IDR_x_test=c()
  for (i in 1:length(y_test)){
    Y=sort(y_test)[i]
    K=c()
    for(j in 1:(length(unique_order_y)-1)){
      
      K=append(K,ifelse(nu!=Inf,integrate(function(u){K_h(Y-u)},
                                          lower = unique_order_y[j],
                                          upper = unique_order_y[j+1],
                                          stop.on.error = FALSE)$value,
                        2^(-1)*(erf((Y-unique_order_y[j])/h)-
                                  erf((Y-unique_order_y[j+1])/h))))}
        #upper=0.5+((Y-unique_order_y[j])/h)*gamma((nu+1)/2)*
        #  hypergeom(0.5,(nu+1)/2,1.5,-((Y-unique_order_y[j])/h)^2/nu)/(sqrt(pi*nu)*gamma(nu/2))
        #lower=0.5+((Y-unique_order_y[j+1])/h)*gamma((nu+1)/2)*
         # hypergeom(0.5,(nu+1)/2,1.5,-((Y-unique_order_y[j+1])/h)^2/nu)/(sqrt(pi*nu)*gamma(nu/2))
        #K=append(K,upper-lower)
      
    
    
    K=append(K,ifelse(nu!=Inf,
                      integrate(function(u){K_h(Y-u)},
                                lower = unique_order_y[length(unique_order_y)],
                              upper = Inf)$value,
                      2^(-1)*(erf((Y-unique_order_y[length(unique_order_y)])/h)+1)))
    #K=append(K,2^(-1)*(erf((Y-unique_order_y[length(unique_order_y)])/h)+1))
    
    smooth_IDR_x_test=append(smooth_IDR_x_test,
                             w%*%K)
    if(progress){pb$tick()}
  }
  #toc()
  return(smooth_IDR_x_test)
  
  
}

smooth_IDR_density=function(y,x,h=(log(length(x))/length(x))^(1/9),nu=2.5,y_test=y,x_test){
  n=length(x)
  x=sort(x)
  fit_idr=idr(y,data.frame(x))
  unique_order_y=sort(unique(y))
  
  if (nu!=Inf){
    K_h=function(u){1/h*(1/(nu*pi)^(0.5))*gamma((nu+1)/2)/gamma(nu/2)*(1+(u/h)^2/nu)^(-(nu+1)/2)}
  }
  else{
    K_h=function(u){1/h*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)}
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
    smooth_IDR_x_test=append(smooth_IDR_x_test,
                             w%*%sapply(unique_order_y,function(u){K_h(Y-u)}
                                        
                             ))}
  return(smooth_IDR_x_test)
  
}



