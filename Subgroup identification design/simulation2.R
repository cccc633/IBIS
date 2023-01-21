source("function1.R")
source("function2.R")
K = 3; J = 4  # K rows and J columns
n <- 10 # sample size
scenario <- as.matrix(read.csv("scenarios.csv", header=F)) ## odd-numbered rows are the means 
                                                           ## even number rows are the standard deviations
n.scenario <- nrow(scenario)/2
ntrial = 10000  ## Number of simulated trial for each scenario

###### Independent Analysis #####
## phi (the statistical threshold) should be calibrated
phi = 2.92; #phi = 3.37
## The following matrixs are used to deposit results
Reject.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Total.result <- matrix(NA, nrow = n.scenario, ncol = 2)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, Independent, K, J, n, scenario[2*i-1,], scenario[2*i,], phi)
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixReject <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixTotal <- matrix(NA, nrow = ntrial, ncol = 2)
  for (trial in 1:ntrial){
    matrixReject[trial,] <- temp[[trial]][[1]]
    matrixTotal[trial,] <- c(temp[[trial]][[2]], temp[[trial]][[3]])
  }
  Reject.result[i,] <- colMeans(matrixReject)   ## Calculate the rejection rate for each subgroup
  Total.result[i,] <- colMeans(matrixTotal)  ## Calculate FWER and Conjunctive power
}
Reject.I <- as.data.frame(Reject.result);         
Reject.I$FWER <- Total.result[,1]; Reject.I$Power <- Total.result[,2]
Reject.I$Scenario <- 1:n.scenario; Reject.I$Method <- "Independent";

###### BHM #####
## phi (the statistical threshold) should be calibrated
phi = 9.90; # phi = 17.0
## The following matrixs are used to deposit results
Reject.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Total.result <- matrix(NA, nrow = n.scenario, ncol = 2)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, BHM, K, J, n, scenario[2*i-1,], scenario[2*i,], phi)
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixReject <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixTotal <- matrix(NA, nrow = ntrial, ncol = 2)
  for (trial in 1:ntrial){
    matrixReject[trial,] <- temp[[trial]][[1]]  
    matrixTotal[trial,] <- c(temp[[trial]][[2]], temp[[trial]][[3]])  
  }
  Reject.result[i,] <- colMeans(matrixReject)  ## Calculate the rejection rate for each subgroup
  Total.result[i,] <- colMeans(matrixTotal)  ## Calculate FWER and Conjunctive power
}
Reject.B <- as.data.frame(Reject.result);         
Reject.B$FWER <- Total.result[,1]; Reject.B$Power <- Total.result[,2]
Reject.B$Scenario <- 1:n.scenario; Reject.B$Method <- "BHM";

###### IBIS #####
## phi (the statistical threshold) should be calibrated
phi = 118.1; # phi = 311.5
## Obtain all possible subgroup divisions
G.futile <- Divide34()
## The following matrixs are used to deposit results
Reject.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Total.result <- matrix(NA, nrow = n.scenario, ncol = 2)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, IBIS, K, J, n, scenario[2*i-1,], scenario[2*i,], G.futile, phi)
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixReject <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixTotal <- matrix(NA, nrow = ntrial, ncol = 2)
  for (trial in 1:ntrial){
    matrixReject[trial,] <- temp[[trial]][[1]]  
    matrixTotal[trial,] <- c(temp[[trial]][[2]], temp[[trial]][[3]])  
  }
  Reject.result[i,] <- colMeans(matrixReject)  ## Calculate the rejection rate for each subgroup
  Total.result[i,] <- colMeans(matrixTotal)  ## Calculate FWER and Conjunctive power
}
Reject.IBIS <- as.data.frame(Reject.result);         
Reject.IBIS$FWER <- Total.result[,1]; Reject.IBIS$Power <- Total.result[,2]
Reject.IBIS$Scenario <- 1:n.scenario; Reject.IBIS$Method <- "IBIS";


###### Freq #####
## phi (the statistical threshold) should be calibrated
phi = 2.29; # phi = 2.66
## Obtain all possible subgroup divisions
G.futile <- Divide34()
## The following matrixs are used to deposit results
Reject.result <- matrix(NA, nrow = n.scenario, ncol = K*J)
Total.result <- matrix(NA, nrow = n.scenario, ncol = 2)
for (i in 1:n.scenario){
  library(parallel)
  ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
  cl <- makeCluster(mc <- getOption("cl.cores", 60))
  simVector <- 1:ntrial
  temp <- clusterApplyLB(cl, simVector, Freq, K, J, n, scenario[2*i-1,], scenario[2*i,], G.futile, phi)
  stopCluster(cl)
  
  ## Read results of simulated trials
  matrixReject <- matrix(NA, nrow = ntrial, ncol = J*K)
  matrixTotal <- matrix(NA, nrow = ntrial, ncol = 2)
  for (trial in 1:ntrial){
    matrixReject[trial,] <- temp[[trial]][[1]]
    matrixTotal[trial,] <- c(temp[[trial]][[2]], temp[[trial]][[3]])
  }
  Reject.result[i,] <- colMeans(matrixReject)  ## Calculate the rejection rate for each subgroup
  Total.result[i,] <- colMeans(matrixTotal)  ## Calculate FWER and Conjunctive power
}
Reject.Freq <- as.data.frame(Reject.result);         
Reject.Freq$FWER <- Total.result[,1]; Reject.Freq$Power <- Total.result[,2]
Reject.Freq$Scenario <- 1:n.scenario; Reject.Freq$Method <- "Freq";


write.csv(rbind(Reject.I, Reject.B, Reject.IBIS, Reject.Freq),"Reject.csv")


