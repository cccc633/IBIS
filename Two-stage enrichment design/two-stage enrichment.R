source("function1.R")
source("function2.R")
K = 3; J = 4  # K rows and J columns
n1 <- 5 # sample size for the 1st stage
n2 <- 5 # sample size for the 2nd stage
scenario <- as.matrix(read.csv("scenarios.csv", header=F)) ## odd-numbered rows are the means 
                                                           ## even number rows are the standard deviations
n.scenario <- nrow(scenario)/2
ntrial = 500  ## Number of simulated trial
Sc = 8  ## Calibrate thresholds for the 8th scenario 
###### IBIS #####
phi.e <- seq(50, 250, by=10)  ## Range for calibrating phi.e
phi.p <- seq(1, 20, by=1)     ## Range for calibrating phi.p
## Obtain all possible subgroup divisions
G.futile <- Divide34()
## The following matrixs are used to deposit results
Result <- matrix(NA, nrow = length(phi.e)*length(phi.p), ncol = 5)
for (i in 1:length(phi.e)){
  for (j in 1:length(phi.p)){
    library(parallel)
    ## Set 60 computer cores to be used in parallel computation. It can be changed based on computing resources
    cl <- makeCluster(mc <- getOption("cl.cores", 60))
    simVector <- 1:ntrial
    temp <- clusterApplyLB(cl, simVector, IBIS, K, J, n1, n2, scenario[2*Sc-1,], scenario[2*Sc,], G.futile, phi.e[i], phi.p[j])
    stopCluster(cl)
    
    ## Read results of simulated trials
    matrixResult <- matrix(NA, nrow = ntrial, ncol = 3)
    for (trial in 1:ntrial){
      matrixResult[trial,] <- c(temp[[trial]])
    }
    Result[length(phi.p)*(i-1)+j, 1:3] <- colMeans(matrixResult) ## Calculate FWER, Conjunctive power and EN
    Result[length(phi.p)*(i-1)+j, 4] <- phi.e[i]
    Result[length(phi.p)*(i-1)+j, 5] <- phi.p[j]
  }
}
Result <- as.data.frame(Result)
names(Result) <- c("FWER", "Power", "EN", "phi.e", "phi.p")

write.csv(Result, "Result.csv")