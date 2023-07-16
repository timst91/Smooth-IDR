library(isodistrreg)
library(progress)
library(tictoc)

library(cubature)


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
  if(time){Time=toc(echo=1)
  
  
  return(list(logS=n/(n-count)*mean(logS),time=unname(Time)))
  }else{
    return(n/(n-count)*mean(logS))
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
    logS[i]=min(c(-log(model),-log(.Machine$double.xmin)))
    
    if(progress){pb$tick()}
  }
  if(progress){pb$terminate()}
  if(time){
    time=toc(echo=1)
    return(list(logS=mean(logS),time=unname(time)))
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
    K_h=function(x,h){dnorm(x,0,h)}
  }else{
    K_h=function(x,h){dt(x/h,nu)/h}
  }
  
  
  if(nu_init!=Inf){
    K_h_2=function(u,h){-1/h^3*(nu_init*pi)^(-0.5)*gamma((nu_init+1)/2)/gamma(nu_init/2)*(nu_init+1)/nu_init*(
      (1+(u/h)^2/nu_init)^(-(nu_init+3)/2)-(nu_init+3)/nu_init*(u/h)^2*(1+(u/h)^2/nu_init)^(-(nu_init+5)/2))}
  }else{
    K_h_2=function(u,h){1/h^3*(2*pi)^(-0.5)*exp(-2*(u/h)^2/2)*(u^2/h^2-1)}
  }
  
  int=smooth_IDR_density_h_opt(y=y,x=x,
                               c=c,nu=nu,
                               y_test = y[1],
                               x_test =x[1],
                               #range_int = Inf,
                               normalize = 1)$integral
  
  
  count=0
  logS=numeric(n)
  
  kappa0=K_h(0,1)
  
  eps_n=(log(n)/n)^1/3
  var=ifelse(nu!=Inf,nu/(nu-2),1)
  
  for (i in 1:n){
    w=cdf[i,1]+numeric(m)
    for (j in 2:m){
      w[j]=cdf[i,j]-cdf[i,j-1]
      
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
      logS[i]=-log(1/int*w_i%*%K)
    }else{
      count=count+1
      logS[i]= ifelse(i>1,1/int*logS[i-1],0)
      }
    
    if(progress){pb$tick()}
  }
  
  if(progress){pb$terminate()}
  count=0
  if(time){
    time=toc(quiet=TRUE)
    return(list(logS=mean(logS),time=unname(time$toc-time$tic)))
  }else{
    return(n/(n-count)*mean(logS))
  }
  
}

CV_c=function(y,x,c,nu,progress=FALSE){
  
  n=length(x)
  
  logS=numeric(n)
  
  if(progress){pb= progress_bar$new(total=length(x))}
  int=smooth_IDR_density_h_opt(y=y,x=x,
                               c=c,nu=nu,
                               y_test = y[1],
                               x_test =x[1],
                              # range_int = Inf,
                               normalize = 1)$integral
  
  for (i in 1:n){
    #model=1/int*
     model=1/int* smooth_IDR_density_h_opt(y=y[-i],x=x[-i],
                                         c=c,nu=nu,
                                         y_test = y[i],
                                         x_test =x[i])$density
    
    logS[i]=min(c(-log(model),-log(.Machine$double.xmin)))
    if(progress){pb$tick()}
  }
  if(progress){pb$terminate()}
  
  
  return(mean(logS))
}

