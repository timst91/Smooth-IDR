library(neuralnet)
library(MASS)
library(progress)


data(Boston)

normalize <- function(x) {
  return((x - mean(x)) / sd(x))
}

target <- normalize(Boston$medv)
data <- data.frame(target, lapply(Boston[, -14], normalize))

nsplits=20
mean_test_logS=0
NU=c(2.01,3,4,5,10,20,Inf)
H=seq(0.01,5,l=10)

grid=expand.grid(NU,H)

pb=progress_bar$new(total=nsplits*nrow(grid))

for (split in 1:nsplits){
  
  train_indices <- sample(1:nrow(data), 0.72 * nrow(data))
  #train_indices <- sample(1:nrow(data), 0.9 * nrow(data))
  
  train_data <- data[train_indices, ]
  
  rem_data=data[-train_indices, ]
  val_ind=sample(1:nrow(rem_data), 0.18 * nrow(data))
  val_data <- rem_data[val_ind, ]
  test_data=rem_data[-val_ind,]
  
  lr=seq(0.01,2,l=10)
  MSE=c()
  nuu=c()
  hh=c()
  for (Lr in lr){
    formula <- as.formula(paste("target ~", paste(names(data)[-1], collapse = " + ")))
    
    model <- neuralnet(formula, data = train_data, hidden = c(50, 50), linear.output = TRUE, 
                       act.fct = "logistic",learningrate = Lr, stepmax = 1e6)
    
    
    predictions <- compute(model, val_data[, -1])$net.result
    MSE=append(MSE,mean((predictions-val_data[, 1])^2))
    
    X=compute(model, train_data[, -1])$net.result
    y=train_data$target
    
    of_h_nu=c()
    
    
    for (i in 1:nrow(grid)){
      
      nu=as.numeric(grid[i,])[1]
      h=as.numeric(grid[i,])[2]
      of=OF_h(y,X,h,nu)$logS
      of_h_nu=append(of_h_nu,
                     of)
      #print(of)
      pb$tick()
    }
    
    nuu=append(nuu,as.numeric(grid[which.min(of_h_nu),])[1])
    hh=append(hh,as.numeric(grid[which.min(of_h_nu),])[2])
  }
  lr=lr[which.min(MSE)]
  h_opt=hh[which.min(MSE)]
  nu_opt=nuu[which.min(MSE)]
  
  
  train_val_ind=c(train_indices,val_ind)
  model <- neuralnet(formula, data = data[train_val_ind,], hidden = c(50, 50), linear.output = TRUE, 
                     act.fct = "logistic", stepmax = 1e6)
  
  
  # Step 7: Make predictions using the trained model
  predictions <- compute(model, data[-train_val_ind, -1])$net.result
  
  predictions <- compute(model, val_data[, -1])$net.result
  
  X=compute(model, train_data[, -1])$net.result
  y=train_data$target
  
  
  
  val_logS=c()
  for (i in 1:nrow(grid)){
    
    nu=as.numeric(grid[i,])[1]
    h=as.numeric(grid[i,])[2]
    pred_idr=c()
    
    for(j in 1:nrow(val_data)){
      pred_idr=append(pred_idr,smooth_IDR_CDF(y,X,h=h,nu=nu,x_test = predictions[j],
                                              y_test = val_data$target[j]))
    }
    
    val_logS=append(val_logS,-sum(log(pred_idr)))
    #print(of)
    pb$tick()
  }
  of_h_nu=val_logS
  nu_opt=as.numeric(grid[which.min(of_h_nu),])[1]
  h_opt=as.numeric(grid[which.min(of_h_nu),])[2]
  
  train_val_ind=c(train_indices,val_ind)
  model <- neuralnet(formula, data = data[train_val_ind,], hidden = c(50, 50), linear.output = TRUE, 
                     act.fct = "logistic", stepmax = 1e6)
  
  
  # Step 7: Make predictions using the trained model
  predictions <- compute(model, data[-train_val_ind, -1])$net.result
  
  X=compute(model, data[train_val_ind, -1])$net.result
  y=data[train_val_ind,1]
  
  
  pred_idr=c()
  
  for(j in 1:length(predictions)){
    pred_idr=append(pred_idr,smooth_IDR_CDF(y,X,h=h_opt,nu=nu_opt,x_test = predictions[j],
                                            y_test = data[-train_val_ind,1][j]))
  }
  
  mean_test_logS=mean_test_logS-1/nsplits*sum(log(pred_idr))
  
  print(sum(log(pred_idr)))

  }

pb$terminate()
