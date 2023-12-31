library(isodistrreg)
library(progress)
library(tictoc)

smooth_IDR_CDF=function(y,x,x_test,
                        y_test=y,
                        h=(log(length(x))/length(x))^(1/9),
                        nu=2.5,
                        progress=FALSE,
                        time=FALSE){
  
  if(time){tic()}
  
  
  sort_ind=match(sort(unique(x)),x)
  
  x=sort(unique(x))
  y=y[sort_ind]
  
  unique_order_y=sort(unique(y))
  n=length(x)
  m=length(unique_order_y)
  
  fit.idr=idr(y,data.frame(x))
  pred.idr=predict(fit.idr,data=data.frame(x=x_test))
  
  #compute cdfs
  points=pred.idr[[1]]$points
  cdfs=pred.idr[[1]]$cdf
  
  cdf=numeric(m)
  
  for(j in 1:m){
    y=unique_order_y[j]
    ind_y=ifelse(points[1]<=y,max(which(points <= y)),0)
    
    cdf[j]=ifelse(ind_y>0,cdfs[ind_y],0)
    
  }
  
  K_int=numeric(m)
  
  y_test=sort(y_test)
  m_test=length(y_test)
  
  smooth_IDR_x_test=numeric(m_test)
  if(progress){pb = progress_bar$new(total=m_test)}
  
  
  for (i in 1:m_test){
    Y=y_test[i]
    
    if(nu== Inf){
      for(j in 1:(m-1)){
        K_int[j]=pnorm(unique_order_y[j+1],Y,h)-pnorm(unique_order_y[j],Y,h)
      }
      K_int[m]=1-pnorm(unique_order_y[m],Y,h)
    }else{
      
      for(j in 1:(m-1)){
        K_int[j]=pt((unique_order_y[j+1]-Y)/h,nu)-pt((unique_order_y[j]-Y)/h,nu)
      }
      
      K_int[m]=1-pt((unique_order_y[m]-Y)/h,nu)
    }
    if(progress){pb$tick()}
    
    smooth_IDR_x_test[i]=cdf%*%K_int
  }
  
  if(time){toc()}
  
  if(progress){pb$terminate()}
  
  return(smooth_IDR_x_test)
}   
    


