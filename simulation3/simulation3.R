source("function3.R")
K = 3; J = 4  # K rows and J columns
n1 <- 5 # sample size
n2 <- 5 # sample size
scenario <- as.matrix(read.csv("scenarios.csv", header=F))
n.scenario <- nrow(scenario)/2
ntrial = 1000
###### IBIS #####
phi.e <- seq(150, 250, by=10)
phi.p <- seq(1, 10, by=1)
G.futile <- Divide34()
Result <- matrix(NA, nrow = length(phi.e)*length(phi.p), ncol = 5)
for (i in 1:length(phi.e)){
  for (j in 1:length(phi.p)){
    library(parallel)
    cl <- makeCluster(mc <- getOption("cl.cores", 60))
    simVector <- 1:ntrial
    temp <- clusterApplyLB(cl, simVector, IBIS, K, J, n1, n2, scenario[1,], scenario[2,], G.futile, phi.e[i], phi.p[j])
    stopCluster(cl)
    
    matrixResult <- matrix(NA, nrow = ntrial, ncol = 3)
    for (trial in 1:ntrial){
      matrixResult[trial,] <- c(temp[[trial]])
    }
    Result[length(phi.p)*(i-1)+j, 1:3] <- colMeans(matrixResult)
    Result[length(phi.p)*(i-1)+j, 4] <- phi.e[i]
    Result[length(phi.p)*(i-1)+j, 5] <- phi.p[j]
  }
}
Result <- as.data.frame(Result)
names(Result) <- c("FWER", "Power", "EN", "phi.e", "phi.p")

write.csv(Result, "Result.csv")