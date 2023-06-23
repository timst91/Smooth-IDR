#install.packages("keras")
#install.packages("tidyverse")

library(tictoc)
library(tensorflow)
library(MASS)
library(keras)
library(tidyverse)
library(progress)
reticulate::use_python("C:\\Apps\\Anaconda3\\python.exe")


data(Boston)
df =as.data.frame(Boston)
data=df
# Normalize the predictor variables
normalized_df = as.data.frame(lapply(df[, -14], scale))

y= df$medv

nsplits=20

NU=c(2.01,3,4,5,10,20,Inf)
c=seq(0.1,4,l=10)

h_init=1.6
nu_init=4

grid=expand.grid(NU,c)

lr=seq(0.001,0.01,l=10)
batch_size=seq(10,300,l=10)
grid2=expand.grid(lr,batch_size)
epochs=40

mean_test_logS=0
total=nsplits*nrow(grid1)*nrow(grid2)
pb=progress_bar$new(total=total)
tic()
step=0
for (split in 1:nsplits){
  ind=sample(1:nrow(data))
  train_ind = ind[1:floor(0.72 * nrow(data))]
  val_ind=ind[(floor(0.72 * nrow(data))+1):floor(0.9 * nrow(data))]
  test_ind=(1:nrow(data))[-c(train_ind,val_ind)]
  
  X_train= normalized_df[train_ind, ]
  y_train=y[train_ind]
  
  X_val=normalized_df[val_ind, ]
  y_val=y[val_ind]
  
  X_train_val=normalized_df[c(train_ind,val_ind), ]
  y_train_val=y[c(train_ind,val_ind)]
  
  X_test=normalized_df[-c(train_ind,val_ind), ]
  y_test=y[-c(train_ind,val_ind)]
  
  
  
  MSE=c()
  nuu=c()
  cc=c()
  for (j in 1:nrow(grid2)){
    Lr=as.numeric(grid2[j,])[1]
    bs=as.numeric(grid2[j,])[2]
    
    model = keras_model_sequential()
    model %>%
      layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
      #layer_dense(units = 50, activation = "relu") %>%
      layer_dense(units = 1, activation = "linear")
    
    
    
    model %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(learning_rate  = Lr)
    )
    
    
    
    # Train model
    model %>% fit(
      x = as.matrix(X_train),
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
    
    of_c_nu=c()
    
    for (i in 1:nrow(grid)){
      
      nu=as.numeric(grid[i,])[1]
      c=as.numeric(grid[i,])[2]
      ofc=OF_c(y_train,fitted,c,nu,
               nu_init = nu_init,h_init = h_init)
      of_c_nu=append(of_c_nu,
                     ofc)
      pb$tick()
      step=step+1
     
      if(step%%100==0){
        print(paste(step," / ", total))}
    }
    
    nuu=append(nuu,as.numeric(grid[which.min(of_c_nu),])[1])
    cc=append(cc,as.numeric(grid[which.min(of_c_nu),])[2])
  }
  lr_opt=as.numeric(grid2[which.min(MSE),])[1]
  bs_opt=as.numeric(grid2[which.min(MSE),])[2]
  c_opt=cc[which.min(MSE)]
  nu_opt=nuu[which.min(MSE)]
  
  model = keras_model_sequential()
  model %>%
    layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
    #layer_dense(units = 50, activation = "relu") %>%
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
  
  
  fitted=as.numeric(model %>% predict(as.matrix(X_train_val),
                                      verbose = 0))
  
  predictions = as.numeric(model %>% predict(as.matrix(X_test),
                                             verbose = 0))
  
  test_logS=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density_h_opt(y_train_val,fitted,c=c_opt,
                                 h_init = h_init,
                                 nu_init=nu_init,
                                 nu=nu_opt,
                                 y_test=u[2],
                                 x_test=u[1]
                                #normalize=TRUE
                                )$density)}))
  
  
  mean_test_logS=mean_test_logS-1/nsplits*test_logS
  time=toc(quiet=TRUE)
  print(test_logS)
  
  
}

pb$terminate()
cat(paste('mean test logS:', mean_test_logS,"| ",
          ' time elapsed:',unname(time$toc)))
