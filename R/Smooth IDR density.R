library(isodistrreg)


smooth_IDR_density=function(y,x,x_test,
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
  
  #compute w_j(x_test)
  
  points=pred.idr[[1]]$points
  cdf=pred.idr[[1]]$cdf
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
  
  if (nu==Inf){
    
    K_h=function(x){dnorm(x,0,h)}
  }
  else{
    K_h=function(x){dt(x/h,nu)/h}
  }
  
 
  
  y_test=sort(y_test)
  m_test=length(y_test)
  
  if(progress){pb = progress_bar$new(total=m_test)}
  
  
  smooth_IDR_x_test=numeric(m_test)
  
  for (j in 1:m_test){
    Y=y_test[j]
    smooth_IDR_x_test[j]=w%*%sapply(unique_order_y,function(u){K_h(Y-u)})
    if(progress){pb$tick()}
    }

  if(time){toc()}
  if(progress){pb$terminate()}
  
  return(smooth_IDR_x_test)
  
}
