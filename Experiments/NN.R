
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

NU=c(2.01,3,4,5,10,20,Inf)
H=seq(0.01,5,l=10)

grid=expand.grid(NU,H)
lr=seq(0.01,2,l=10)

#pb=progress_bar$new(total=nsplits*nrow(grid))
mean_test_logS=0

for (split in 1:nsplits){
  ind=sample(1:nrow(data))
  train_ind <- ind[1:floor(0.72 * nrow(data))]
  val_ind=ind[(floor(0.72 * nrow(data))+1):floor(0.9 * nrow(data))]
  test_ind=(1:nrow(data))[-c(train_ind,val_ind)]
  
  train_data <- data[train_ind, ]
  val_data=data[val_ind,]
  test_data=data[test_ind,]
 
  
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
    pb$tick()
    
    for (i in 1:nrow(grid)){
      
      nu=as.numeric(grid[i,])[1]
      h=as.numeric(grid[i,])[2]
      of=OF_h(y,X,h,nu)$logS
      of_h_nu=append(of_h_nu,
                     of)
      #print(of)
      #pb$tick()
    }
    
    nuu=append(nuu,as.numeric(grid[which.min(of_h_nu),])[1])
    hh=append(hh,as.numeric(grid[which.min(of_h_nu),])[2])
  }
  lr=lr[which.min(MSE)]
  h_opt=hh[which.min(MSE)]
  nu_opt=nuu[which.min(MSE)]
  
  
  train_val_ind=c(train_ind,val_ind)
  model <- neuralnet(formula, data = data[train_val_ind,], hidden = c(50, 50), linear.output = TRUE, 
                     act.fct = "logistic", stepmax = 1e6,
                     learningrate = lr)
  
  
  # Step 7: Make predictions using the trained model
  X=as.numeric(model$net.result[[1]])
  y=data[train_val_ind,1] 
  predictions <- compute(model, data[test_ind, -1])$net.result
  
  
  X_test=as.numeric(predictions)
  y_test=test_data$target
  test_logS=mean(apply(cbind(X_test,y_test),1, function(u){
    log(smooth_IDR_density(y,X,h=h_opt,nu=nu_opt,y_test=u[2],x_test=u[1]))}))
  
  
  mean_test_logS=mean_test_logS-1/nsplits*test_logS
  print(test_logS)

  }

#pb$terminate()
print(mean_test_logS)


j=10
smooth_IDR_density(y,X,h=h_opt,nu=nu_opt,y_test=yy,x_test=X_test[j])
c(X_test[j],y_test[j])

plot(yy,smooth_IDR_density(y,X,x_test=X_test[j],y_test = yy),type='l')
abline(v=y_test[j])
