# Install the required packages
#install.packages("keras")
#install.packages("tidyverse")
library(tensorflow)
library(MASS)
library(keras)
library(tidyverse)
reticulate::use_python("C:\\Apps\\Anaconda3\\envs\\myenv\\python.exe")
# Load the Boston Housing dataset from the MASS package
data(Boston)
df <- as.data.frame(Boston)

# Normalize the predictor variables
normalized_df <- as.data.frame(lapply(df[, -14], scale))

y= df$medv

nsplits=20

NU=c(2.01,3,4,5,10,20,Inf)
H=seq(0.01,5,l=10)

grid=expand.grid(NU,H)

lr=seq(0.001,0.01,l=10)
batch_size=seq(50,300,l=10)
grid2=expand.grid(lr,batch_size)
epochs=40

mean_test_logS=0

pb=progress_bar$new(total=nsplits*nrow(grid)*nrow(grid2))

for (split in 1:nsplits){
  ind=sample(1:nrow(data))
  train_ind <- ind[1:floor(0.72 * nrow(data))]
  val_ind=ind[(floor(0.72 * nrow(data))+1):floor(0.9 * nrow(data))]
  test_ind=(1:nrow(data))[-c(train_ind,val_ind)]
  
  X_train<- normalized_df[train_ind, ]
  y_train=y[train_ind]
  
  X_val=normalized_df[val_ind, ]
  y_val=y[val_ind]
  
  X_train_val=normalized_df[c(train_ind,val_ind), ]
  y_train_val=y[c(train_ind,val_ind)]
  
  X_test=normalized_df[-c(train_ind,val_ind), ]
  y_test=y[-c(train_ind,val_ind)]
  

  
  MSE=c()
  nuu=c()
  hh=c()
  for (j in 1:nrow(grid2)){
    Lr=as.numeric(grid2[j,])[1]
    bs=as.numeric(grid2[j,])[2]
    
    model <- keras_model_sequential()
    model %>%
      layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
      layer_dense(units = 50, activation = "relu") %>%
      layer_dense(units = 1, activation = "linear")
    
    # Compile the model
   
    model %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_sgd(learning_rate  = Lr)
    )
    
    
    
    # Train the model
    model %>% fit(
      x = as.matrix(X_train),
      y = y_train,
      batch_size = bs,
      epochs = epochs,
      verbose = 0,
      validation_split = 0
    )
    
    # Evaluate the model on the test set
    MSE=append(MSE,unname(model %>% evaluate(
      x = as.matrix(X_val),
      y = y_val,
      verbose = 0
    )))
    fitted=as.numeric( model %>% predict(as.matrix(X_train),
                                         verbose = 0))
    
    of_h_nu=c()
    #pb$tick()
    
    for (i in 1:nrow(grid)){
      
      nu=as.numeric(grid[i,])[1]
      h=as.numeric(grid[i,])[2]
      of=OF_h(y_train,fitted,h,nu)$logS
      of_h_nu=append(of_h_nu,
                     of)
      #print(of)
      pb$tick()
    }
    
    nuu=append(nuu,as.numeric(grid[which.min(of_h_nu),])[1])
    hh=append(hh,as.numeric(grid[which.min(of_h_nu),])[2])
  }
  lr_opt=as.numeric(grid2[which.min(MSE),])[1]
  bs_opt=as.numeric(grid2[which.min(MSE),])[2]
  h_opt=hh[which.min(MSE)]
  nu_opt=nuu[which.min(MSE)]
  
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 1, activation = "linear")
  
  # Compile the model
  
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_sgd(learning_rate  = lr_opt)
  )
  
  # Set the batch size and number of epochs
  
  
  
  # Train the model
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
 
  predictions <- as.numeric(model %>% predict(as.matrix(X_test),
                                              verbose = 0))
  
  test_logS=mean(apply(cbind(predictions,y_test),1, function(u){
    log(smooth_IDR_density(y_train_val,fitted,h=h_opt,nu=nu_opt,y_test=u[2],x_test=u[1]))}))
  
  
  mean_test_logS=mean_test_logS-1/nsplits*test_logS
  print(test_logS)
  
}

pb$terminate()
print(mean_test_logS)

