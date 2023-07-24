# SUMMARY OF PROCEDURE
#
#
# This code compares the performance between smooth IDR with global and local
# bandwidth selection based on algorithm 1 from 
# https://arxiv.org/pdf/2212.08376v1.pdf. The procedure is as follows:
# For (lr,bs) in (learning rates, batch sizes){
#   train neural network with parameters (lr,bs) 
#   select optimal smooth IDR parameters (nu,h,h_init,c,nu_init) in 4 different 
#   ways:
#    1. For smooth IDR with global bandwidth choose (h,nu) by minimizing 
#       OF(y_val,prediction(x_val),h,nu)
#    2. For smooth IDR with global bandwidth choose (h,nu) by minimizing 
#       -sum(log(smooth_IDR_density_global_h(y_train,prediction(x_train),
#                                           y_eval=y_val,
#                                           x_eval=prediction(x_val),
#                                           h,nu))) 
#    3. For smooth IDR with local bandwidth choose (h_init,nu_init) obtained
#       by optimal (h,nu) from 1. and select (c,nu) by minimizing 
#
#       -sum(log(smooth_IDR_density_local_h(y_train,prediction(x_train),
#                                           y_eval=y_val,
#                                           x_eval=prediction(x_val),
#                                           h_init=h_opt1,nu_init=nu_opt1,
#                                           c,nu))) 
#    4. For smooth IDR with local bandwidth choose (c,nu,h_init,nu_init=nu)
#       by minimizing 
#
#       -sum(log(smooth_IDR_density_local_h(y_train,prediction(x_train),
#                                           y_eval=y_val,
#                                           x_eval=prediction(x_val),
#                                           h_init,nu_init=nu,
#                                           c,nu))) 
#}
# Choose optimal (lr,bs) based on minimal MSE and the according optimal 
# (h,nu,c,nu_init,h_init). For evaluation of test set performance compute
#
#     -sum(log(smooth_IDR(y[training+validation set],
#                         prediction(x[training+validation set]),
#                         y_eval=y_test,x_eval=prediction(x_test))
#                    
# with the parameters chosen with each of the steps 1.-4.
#
# Also, the chosen parameters and normalizing integrals of the local density 
# estimates are being recorded.

###############################################################################


# IMPORTANT INFO REGARDING IMPLEMENTATION
# 
# The code requires the use of a python environment through reticulate::use_python().
# The environment needs to have tensorflow installed. 
# Defining the neural network depends on the versions of python and R; if 
#
# model = keras_model_sequential()
# model %>%
#  layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
#  layer_dense(units = 1, activation = "linear")
#
# gives an error, you need to use tensor notation like: 
#
# model <- tf$keras$Sequential()
# model$add(tf$keras$layers$Dense(units = 50, activation = "relu", 
#                                input_dim = as.integer(13)))
# model$add(tf$keras$layers$Dense(units =50, activation = "relu"))
# model$add(tf$keras$layers$Dense(units = 1, activation = "linear"))
#
# model$compile(
# loss = "mean_squared_error",
# optimizer = tf$keras$optimizers$Adam(learning_rate = Lr)
# )

# Train the model
# model$fit(
#  x = np$array(as.matrix(unname(X_train))),
#  y = np$array(y_train),
# batch_size = as.integer(bs),
#  epochs = as.integer(epochs),
#  verbose = 0,
#  validation_split = 0
# )
#
#
#
# MSE = append(MSE, unname(model$evaluate(
#  x = np$array(as.matrix(unname(X_val))),
#  y = np$array(y_val),
#  verbose = 0
# )))

# fitted <- as.numeric(model$ predict(
#  tf$constant(as.matrix(X_train)), verbose = 0))

###############################################################################





#install.packages("devtools")
library(devtools)
#loads the neccessary functions
source_url("https://raw.githubusercontent.com/timst91/Smooth-IDR/main/R/all.R")



#install.packages("keras")
#install.packages("tidyverse")
#install.packages("cubature")
#install.packages("tensorflow")
#install.packages("tictoc")
#install.packages("progress")
#install.packages("isodistrreg")
#install.packages("reticulate")


library(cubature)
library(tictoc)
library(tensorflow)
library(MASS)
library(keras)
library(tidyverse)
library(progress)
library(isodistrreg)
library(reticulate)

use_python("/opt/anaconda3/envs/myenv2/bin/python3")

data(Boston)
df =as.data.frame(Boston)
data=df
# Normalize the predictor variables
normalized_df = as.data.frame(lapply(df[, -14], scale))

y= df$medv

nsplits=7

NU=c(2.01,3,4,5,10,20,Inf)

C=seq(0.01,4,l=10)

H=seq(0.01,3,l=10)

H_init=seq(0.005,2,l=7)

grid=expand.grid(NU,C)

grid3=expand.grid(NU,H)

#grid4=expand.grid(NU,NU,H_init,C)
grid4=expand.grid(NU,H_init,C)

#4000 epochs:
lr=seq(1e-6,1e-4,l=3)

#400 epochs:
#lr=seq(1e-5,1e-2,l=3)

#40 epochs:
#lr=seq(1e-3,1e-1,l=3)


batch_size=c(1,8,16,32)



grid2=expand.grid(lr,batch_size)

epochs=4000


mean_test_logS_of=numeric(nsplits)


total=nsplits*nrow(grid2)*(nrow(grid3))

progress=1
pb=progress_bar$new(total=total)
tic() 

step=0
params=matrix(nrow=nsplits,ncol=4)

for (split in 2:nsplits){
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
  
  hh=c()
  
  for (j in 1:nrow(grid2)){
    Lr=as.numeric(grid2[j,])[1]
    bs=as.numeric(grid2[j,])[2]
    
    model = keras_model_sequential()
    model %>%
      layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
      layer_dense(units = 50, activation = "relu") %>%
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
      verbose = 1,
      validation_split = 0
    )
    
    
    MSE=append(MSE,unname(model %>% evaluate(
      x = as.matrix(X_val),
      y = y_val,
      verbose = 0
    )))
    
    
    pred_val=as.numeric( model %>% predict(as.matrix(X_train), verbose = 0))
    
    print(1)
    
    of_h_nu=numeric(nrow(grid3))
    
    grid_search=numeric(nrow(grid3))
    
    for (i in 1:nrow(grid3)){
      
      nu=as.numeric(grid3[i,])[1]
      h=as.numeric(grid3[i,])[2]
      of=OF_h(y_train,pred_val,h,nu,progress = 1)
      
      of_h_nu[i]=of
      if(progress){pb$tick()}
      step=step+1
      if(step%%200==0){print(paste(step," / ", total))}
      
    }
    print(2)
    nu_opt_init=as.numeric(grid3[which.min(of_h_nu),])[1]
    
    
    h_opt_init=as.numeric(grid3[which.min(of_h_nu),])[2]
    nuu=append(nuu,nu_opt_init)
    hh=append(hh,h_opt_init)
  }
  
  lr_opt=as.numeric(grid2[which.min(MSE),])[1]
  bs_opt=as.numeric(grid2[which.min(MSE),])[2]
  
  h_opt=hh[which.min(MSE)]
  nu_opt=nuu[which.min(MSE)]
  
  params[split,]=c(lr_opt,bs_opt,h_opt,nu_opt)
  
  model = keras_model_sequential()
  model %>%
    layer_dense(units = 50, activation = "relu", input_shape = 13) %>%
     layer_dense(units = 50, activation = "relu") %>%
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
  
  
  
  
  
  mean_test_logS_of[split]=test_logS_of
  
  
  
  
  print(-c(test_logS_of))
  
  
}

time=toc()
pb$terminate()
cat(paste('mean test logS of global:',-mean( mean_test_logS_of)))

## integrals: 

normalizing_integrals

## parameters:

params

colnames(params)= c("lr_opt","bs_opt","c_opt","nu_opt_local","nu_opt_gs", "nu_opt_of",
                    "h_opt_gs", "h_opt_of","nu2_opt", "h_opt_init","c_opt_2")

res=-cbind(mean_test_logS_local, mean_test_logS_local_norm,
           mean_test_logS_local2, mean_test_logS_local_norm2,
           mean_test_logS_of,mean_test_logS_grid_search)
colnames(res)=c("local bw 1","local bw 1 normalized",
                "local bw 2","local bw 2 normalized",
                "global bw of", "global bw gs")




# Specify the file path for the export
file_path = "/Users/ruby.mars/Desktop/Params_experiment/params400_1.txt"

# Export the matrix as a .txt file
write.table(params, file = file_path, sep = "\t", row.names = FALSE,
            quote = FALSE)

# Specify the file path for the export
file_path = "/Users/ruby.mars/Desktop/Params_experiment/logS_scores_400_1.txt"

# Export the matrix as a .txt file
write.table(res, file = file_path, sep = "\t", row.names = FALSE,
            quote = FALSE)

data=read.table(file_path,header=TRUE,sep="\t")

data

apply(data,2,mean)
