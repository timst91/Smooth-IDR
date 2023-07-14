#install.packages("pracma")
library(pracma)

#install.packages("keras")
#install.packages("tidyverse")

library(tictoc)
library(tensorflow)
library(MASS)
library(keras)
library(tidyverse)
library(progress)
library(isodistrreg)
library(reticulate)

reticulate::use_python("/opt/anaconda3/envs/myenv2/bin/python3")

data(Boston)
df =as.data.frame(Boston)
data=df
# Normalize the predictor variables
normalized_df = as.data.frame(lapply(df[, -14], scale))

y= df$medv

nsplits=1

NU=c(2.01,3,4,5,10,20,Inf)

C=seq(0.01,4,l=10)

H=seq(0.01,3,l=10)

H_init=seq(0.01,2,l=7)

grid=expand.grid(NU,C)

grid3=expand.grid(NU,H)

#grid4=expand.grid(NU,NU,H_init,C)
grid4=expand.grid(NU,H_init,C)

#4000 epochs:
#lr=seq(1e-6,1e-4,l=3)

#400 epochs:
#lr=seq(1e-5,1e-2,l=3)

#40 epochs:
lr=seq(1e-3,1e-1,l=3)


batch_size=c(1,8,16,32)


grid2=expand.grid(lr,batch_size)

epochs=40

mean_test_logS_local=numeric(nsplits)
mean_test_logS_local_norm=numeric(nsplits)

mean_test_logS_local2=numeric(nsplits)
mean_test_logS_local_norm2=numeric(nsplits)

mean_test_logS_of=numeric(nsplits)
mean_test_logS_grid_search=numeric(nsplits)


total=nsplits*nrow(grid2)*(nrow(grid)+nrow(grid3)+nrow(grid4))

progress=1
pb=progress_bar$new(total=total)
tic() 

step=0

params=matrix(nrow=nsplits,ncol=11)


for (split in 1:nsplits){
  ind=sample(1:nrow(data))
  train_ind = ind[1:floor(0.72 * nrow(data))]
  val_ind=ind[(floor(0.72 * nrow(data))+1):floor(0.9 * nrow(data))]
  test_ind=(1:nrow(data))[-c(train_ind,val_ind)]
  
  X_train= normalized_df[train_ind,]
  y_train=y[train_ind]
  
  X_val=normalized_df[val_ind,]
  y_val=y[val_ind]
  
  X_train_val=normalized_df[c(train_ind,val_ind), ]
  y_train_val=y[c(train_ind,val_ind)]
  
  X_test=normalized_df[-c(train_ind,val_ind), ]
  y_test=y[-c(train_ind,val_ind)]
  
  
  
  MSE=c()
  nuu=c()
  
  nu2=c()
  
  nuu2=c()
  hh=c()
  
  hh_gs=c()
  nuu_gs=c()
  
  cc=c()
  
  hh_init=c()
  cc2=c()
  
  for (j in 1:nrow(grid2)){
    Lr=as.numeric(grid2[j,])[1]
    bs=as.numeric(grid2[j,])[2]
    
    model = keras_model_sequential()
    model %>%
      layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
      # layer_dense(units = 50, activation = "relu") %>%
      layer_dense(units = 1, activation = "linear")
    
    
    
    model %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(learning_rate  = Lr)
    )
    
    
    # Train model
    model %>% fit(
      x =  as.matrix(X_train),
      y = y_train,
      batch_size = bs,
      epochs = epochs,
      verbose = 0,
      validation_split = 0
    )
    
    
    MSE=append(MSE,unname(model %>% evaluate(
      x = as.matrix(X_val),
      y = y_val,
      verbose = 0
    )))
    
    fitted=as.numeric( model %>% predict(as.matrix(X_train),
                                         verbose = 0))
    
    pred_val=as.numeric( model %>% predict(as.matrix(X_val), verbose = 0))
    
    print(1)
    
    of_c_nu=numeric(nrow(grid))
    
    of_c_h_init=numeric(nrow(grid4))
    
    of_h_nu=numeric(nrow(grid3))
    
    grid_search=numeric(nrow(grid3))
    
    for (i in 1:nrow(grid3)){
      
      nu=as.numeric(grid3[i,])[1]
      h=as.numeric(grid3[i,])[2]
      of=OF_h(y_val,pred_val,h,nu)
      
      of_h_nu[i]=of
      
      gs=-sum(apply(cbind(pred_val,y_val),1, function(u){
        log(smooth_IDR_density(y_train,fitted,h=h,nu=nu,
                               y_test=u[2],x_test=u[1]))}))
      
      ### gs should perform better as it optimizes LogS
      ### of density((y,x)_train,(y,x)_test)
      
      grid_search[i]=gs
      
      if(progress){pb$tick()}
      step=step+1
      if(step%%200==0){print(paste(step," / ", total))}
      
    }
    print(2)
    nu_opt_init=as.numeric(grid3[which.min(of_h_nu),])[1]
    
    
    h_opt_init=as.numeric(grid3[which.min(of_h_nu),])[2]
    
    hh=append(hh,h_opt_init)
    nuu2=append(nuu2,nu_opt_init)
    
    nuu_gs=append(nuu_gs,as.numeric(grid3[which.min(grid_search),])[1])
    hh_gs=append(hh_gs,as.numeric(grid3[which.min(grid_search),])[2])
    
    
    for (i in 1:nrow(grid)){
      
      nu=as.numeric(grid[i,])[1]
      c=as.numeric(grid[i,])[2]
      
      val_logS=-sum(apply(cbind(pred_val,y_val),1, function(u){
        log(smooth_IDR_density_h_opt(y_train,fitted,y_test=u[2],x_test=u[1],
                                     nu=nu,c=c,
                                     nu_init = nu_opt_init,
                                     h_init = h_opt_init)$density)}))
      
      of_c_nu[i]=val_logS
      
      if(progress){pb$tick()}
      step=step+1
      if(step%%200==0){print(paste(step," / ", total))}
      
    }
    print(3)
    
    nuu=append(nuu,as.numeric(grid[which.min(of_c_nu),])[1])
    cc=append(cc,as.numeric(grid[which.min(of_c_nu),])[2])
    
    for (i in 1:nrow(grid4)){
      
      nu=as.numeric(grid4[i,])[1]
    
      
      h_init=as.numeric(grid4[i,])[2]
      c=as.numeric(grid4[i,])[3]
      
      val_logS=-sum(apply(cbind(pred_val,y_val),1, function(u){
        log(smooth_IDR_density_h_opt(y_train,fitted,
                                     y_test=u[2],
                                     x_test=u[1],
                                     nu=nu,
                                     c=c,
                                     h_init = h_init)$density)}))
      
      of_c_h_init[i]=val_logS
      
      if(progress){pb$tick()}
      step=step+1
      if(step%%200==0){print(paste(step," / ", total))}
      
    }
    print(4)
    nu2 =append(nu2,as.numeric(grid4[which.min(of_c_h_init),])[1])
    
    
    hh_init =append(hh_init,as.numeric(grid4[which.min(of_c_h_init),])[2])
    cc2=append(cc2,as.numeric(grid4[which.min(of_c_h_init),])[3])
    
  }
  
  lr_opt=as.numeric(grid2[which.min(MSE),])[1]
  bs_opt=as.numeric(grid2[which.min(MSE),])[2]
 
  c_opt=cc[which.min(MSE)]
  nu_opt_local=nuu[which.min(MSE)]
  
  h_opt=hh[which.min(MSE)]
  nu_opt=nuu2[which.min(MSE)]
  
  h_opt_gs=hh_gs[which.min(MSE)]
  nu_opt_gs=nuu_gs[which.min(MSE)]
  
  nu2_opt=nu2[which.min(MSE)]

  h_init_opt=hh_init[which.min(MSE)]
  c_opt_2=cc2[which.min(MSE)]
  
  
 # print(paste('opt lr:', lr_opt,"| ",
  #            ' bs:',bs_opt,"| ",
  #            ' c:', c_opt, "| ",
  #            ' nu local:', nu_opt_local, "| ",
  #            ' nu gs:', nu_opt_gs, "| ",
#              ' nu of:', nu_opt, "| ",
 #             ' h gs:', h_opt_gs, "| ",
  #            ' h of:', h_opt, "| ",
#              '2nd nu:',nu2_opt,"| ",
  #            '2nd nu_init:',nu_init2_opt,"| ",
 #             
   #           'opt h_init:',h_init_opt, "| ",
    #          '2nd c:',c_opt_2,"| ",
  #))
  
  params[split,]=c(lr_opt,bs_opt,c_opt,nu_opt_local,nu_opt_gs, nu_opt,
                   h_opt_gs, h_opt,nu2_opt, h_opt_init,c_opt_2)
  
  
  model = keras_model_sequential()
  model %>%
    layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
    # layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 1, activation = "linear")
  
  
  
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_adam(learning_rate  = lr_opt)
  )
  
  model %>% fit(
    x = as.matrix(X_train_val),
    y = y_train_val,
    batch_size = bs_opt,
    epochs = epochs,
    validation_split = 0,
    verbose = 0
  )
  
  print(5)
  fitted=as.numeric(model %>% predict(as.matrix(X_train_val),
                                      verbose = 0))
  
  predictions = as.numeric(model %>% predict(as.matrix(X_test),
                                             verbose = 0))
  
  
  test_logS_of=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density(y_train_val,fitted,h=h_opt,nu=nu_opt,
                           y_test=u[2],x_test=u[1]))}))
  
  test_logS_grid_search=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density(y_train_val,fitted,h=h_opt_gs,nu=nu_opt_gs,
                           y_test=u[2],x_test=u[1]))}))
  
  
  test_logS_local=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density_h_opt(y_train_val,fitted,
                                 c=c_opt,
                                 h_init = h_opt,
                                 nu_init=nu_opt,
                                 nu=nu_opt_local,
                                 y_test=u[2],
                                 x_test=u[1]
                                 #normalize=TRUE
    )$density)}))
  
  test_logS_local_norm=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density_h_opt(y_train_val,fitted,
                                 c=c_opt,
                                 h_init = h_opt,
                                 nu_init=nu_opt,
                                 nu=nu_opt_local,
                                 y_test=u[2],
                                 x_test=u[1],
                                 normalize=TRUE
    )$density)}))
  
  test_logS_local2=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density_h_opt(y_train_val,fitted,
                                 c=c_opt_2,
                                 h_init = h_init_opt,
                                 nu=nu2_opt,
                                 y_test=u[2],
                                 x_test=u[1]
                                 #normalize=TRUE
    )$density)}))
  
  test_logS_local_norm2=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density_h_opt(y_train_val,fitted,
                                 c=c_opt_2,
                                 h_init = h_init_opt,
                                 nu=nu2_opt,
                                 y_test=u[2],
                                 x_test=u[1],
                                 normalize=TRUE
    )$density)}))
  
  
  
  
  mean_test_logS_local[split]=test_logS_local
  
  mean_test_logS_local_norm[split] =test_logS_local_norm
  
  mean_test_logS_of[split]=test_logS_of
  
  mean_test_logS_grid_search[split]=test_logS_grid_search
  
  mean_test_logS_local2[split]=test_logS_local2
  
  mean_test_logS_local_norm2[split]=test_logS_local_norm2
  
  
  
  print(-c(test_logS_local
           ,test_logS_local_norm,
           test_logS_of,
           test_logS_grid_search,
           test_logS_local2,
           test_logS_local_norm2))
  
  
}

time=toc(echo=0)
pb$terminate()
cat(paste('mean test logS local bw:', -mean(mean_test_logS_local),"| ",
          'mean test logS local bw normal:',-mean( mean_test_logS_local_norm),"| ",
          'mean test logS local bw 2:', -mean(mean_test_logS_local2),"| ",
          'mean test logS local bw normal 2:', -mean(mean_test_logS_local_norm2),"| ",
          
          'mean test logS of global:',-mean( mean_test_logS_of),"| ",
          'mean test logS grid search:', -mean(mean_test_logS_grid_search),"| ",
          
          ' time elapsed:',unname(time)))




