source("function1.R")
source("function2.R")
K = 3; J = 4  # K rows and J columns
n <- 10 # sample size
scenario <- as.matrix(read.csv("scenarios.csv", header=F)) ## odd-numbered rows are the means 
                                                           ## even number rows are the standard deviations
n.scenario <- nrow(scenario)/2
ntrial = 10000  ## Number of simulated trial for each scenario

###### Independent Analysis #####
## These three matrixs are used to deposit results
Bias.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
MSE.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Width.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))  
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, Independent, K, J, n, scenario[2*i-1,], scenario[2*i,])
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixBias <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixMSE <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixWidth <- matrix(NA, nrow = ntrial, ncol = J*K)
  for (trial in 1:ntrial){
    matrixBias[trial,] <- temp[[trial]][[1]]
    matrixMSE[trial,] <- temp[[trial]][[2]]
    matrixWidth[trial,] <- temp[[trial]][[3]]
  }
  ## Calculate the average bias, average MSE, average 95%CI width
  Bias.result[i,] <- colMeans(matrixBias)  
  MSE.result[i,] <- colMeans(matrixMSE)  
  Width.result[i,] <- colMeans(matrixWidth)
}
Bias.I <- as.data.frame(Bias.result);         
MSE.I <- as.data.frame(MSE.result);           
Width.I <- as.data.frame(Width.result);       
Bias.I$Scenario <- 1:n.scenario;     Bias.I$Method <- "Independent";
MSE.I$Scenario <- 1:n.scenario;      MSE.I$Method <- "Independent";
Width.I$Scenario <- 1:n.scenario;    Width.I$Method <- "Independent"

###### BHM #####
## These three matrixs are used to deposit results
Bias.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
MSE.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Width.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, BHM, K, J, n, scenario[2*i-1,], scenario[2*i,])
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixBias <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixMSE <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixWidth <- matrix(NA, nrow = ntrial, ncol = J*K)
  for (trial in 1:ntrial){
    matrixBias[trial,] <- temp[[trial]][[1]]
    matrixMSE[trial,] <- temp[[trial]][[2]]
    matrixWidth[trial,] <- temp[[trial]][[3]]
  }
  ## Calculate the average bias, average MSE, average 95%CI width
  Bias.result[i,] <- colMeans(matrixBias)  
  MSE.result[i,] <- colMeans(matrixMSE)  
  Width.result[i,] <- colMeans(matrixWidth)
}
Bias.B <- as.data.frame(Bias.result);         
MSE.B <- as.data.frame(MSE.result);           
Width.B <- as.data.frame(Width.result);       
Bias.B$Scenario <- 1:n.scenario;     Bias.B$Method <- "BHM";
MSE.B$Scenario <- 1:n.scenario;      MSE.B$Method <- "BHM";
Width.B$Scenario <- 1:n.scenario;    Width.B$Method <- "BHM"

###### IBIS #####
G.futile <- Divide34()  ## Obtain all possible subgroup divisions
## These three matrixs are used to deposit results
Bias.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
MSE.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Width.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, IBIS, K, J, n, scenario[2*i-1,], scenario[2*i,], G.futile)
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixBias <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixMSE <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixWidth <- matrix(NA, nrow = ntrial, ncol = J*K)
  for (trial in 1:ntrial){
    matrixBias[trial,] <- temp[[trial]][[1]]
    matrixMSE[trial,] <- temp[[trial]][[2]]
    matrixWidth[trial,] <- temp[[trial]][[3]]
  }
  ## Calculate the average bias, average MSE, average 95%CI width
  Bias.result[i,] <- colMeans(matrixBias)  
  MSE.result[i,] <- colMeans(matrixMSE)  
  Width.result[i,] <- colMeans(matrixWidth)
}
Bias.IBIS <- as.data.frame(Bias.result);         
MSE.IBIS <- as.data.frame(MSE.result);           
Width.IBIS <- as.data.frame(Width.result);       
Bias.IBIS$Scenario <- 1:n.scenario;     Bias.IBIS$Method <- "IBIS";
MSE.IBIS$Scenario <- 1:n.scenario;      MSE.IBIS$Method <- "IBIS";
Width.IBIS$Scenario <- 1:n.scenario;    Width.IBIS$Method <- "IBIS"

## Arrange results and save results
Bias <- rbind(Bias.I, Bias.B, Bias.IBIS); Bias <- Bias[order(Bias$Scenario),]
write.csv(Bias, "Bias.csv")
MSE <- rbind(MSE.I, MSE.B, MSE.IBIS); MSE <- MSE[order(MSE$Scenario),]
write.csv(MSE, "MSE.csv")
Width <- rbind(Width.I, Width.B, Width.IBIS); Width <- Width[order(Width$Scenario),]
write.csv(Width, "Width.csv")




