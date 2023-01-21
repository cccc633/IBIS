IBIS <- function(trial, K, J, n1, n2, mu, sigma, G.futile, phi.e, phi.p){
  ## This function is used to generate two-stage simulated trial with IBIS analysis ##
  ## Biomarker1 has K levels; Biomarker2 has J levels
  ## n1, sample size for the 1st stage; n2, sample size for the 2nd stage
  ## mu, mean treatment effects for each subgroup
  ## sigma, the standard deviation of treatment effect for each subgroup
  ## G.futile, all possible low-efficacy subsets (obtained by “Divide34” function in “function.R") 
  ## phi.e and phi.p, the statistical threshold
  source("function1.R")
  ##### stage 1 #####
  set.seed(trial)
  data <- matrix(NA, nrow = n1, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n1, mean = mu[i], sd = sigma[i])
  }
  ## Calculate JSD for all possible subgroup divisions
  dis <- list()
  for (i in 1:length(G.futile)){
    matrixfutile <- G.futile[[i]]
    if (i == 1 | i == J*K-1){
      ## When the size of G.futile == 1, there is only one possible situation and we just need to calculate one JSD
      temp.futile <- matrixfutile
      temp.eff <- setdiff(seq(1,J*K), temp.futile)
      dis[[i]] <- calculate.Dnew(data, temp.futile, temp.eff)
    } else {
      result.dis <- c()
      for (j in 1:nrow(matrixfutile)){
        ## When the size of G.futile > 1, there is more than one possible situations and we need to calculate nrow(matrixfutile) JSDs
        temp.futile <- matrixfutile[j,]
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        result.dis[j] <- calculate.Dnew(data, temp.futile, temp.eff)
      }
      dis[[i]] <- result.dis
    }
  }
  ## Target the optimal subgroup division with JSD.max
  dis <- as.numeric(unlist(dis))
  dis.max <- which(dis == max(dis))[1] ## Obtain the serial number of the optimal subgroup division
  for (i in 1:length(G.futile)){ ## Find the high-efficacy(cluster.eff) and low-efficacy(cluster.futile) subset by recursion
    matrixfutile <- G.futile[[i]]
    if (i == 1 | i == J*K-1){
      temp.futile <- matrixfutile
      temp.eff <- setdiff(seq(1,J*K), temp.futile)
      dis.max <- dis.max - 1
      if (dis.max == 0){cluster.futile <- temp.futile; cluster.eff <- temp.eff; break}
    }
    else {
      for (j in 1:nrow(matrixfutile)){
        temp.futile <- matrixfutile[j,]
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        dis.max <- dis.max - 1
        if (dis.max == 0){cluster.futile <- temp.futile; cluster.eff <- temp.eff; break}
      }
    }
  }
  ## inference
  result <- c()
  result0 <- inference.1st(data, cluster.futile, cluster.eff, phi.e, phi.p) ## Call the inference.1st function in "function1.R"
  result[cluster.eff] <- result0[[1]]
  result[cluster.futile] <- result0[[2]]
  ## Make the inferences conform to the monotonicity assumption
  result0 <- matrix(result, nrow = K, byrow = TRUE)
  mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
  if (sum(result) > 0){  
    half <- which(result == 0.5)
    for (i in half){
      place <- which(mm == i, arr.ind = TRUE) ## Target the subgroup with promising conclusion
      result0[place[1]:K, place[2]:J] <- 0.5
    }
    one <- which(result == 1)
    for (i in one){
      place <- which(mm == i, arr.ind = TRUE) ## Target the subgroup with effective conclusion
      result0[place[1]:K, place[2]:J] <- 1
    }
    result <- as.numeric(t(result0))
  }
  N.actual <- rep(n1, J*K)  ## actually enrolled sample size
  ##### stage 2 #####
  if (sum(result == 0.5) > 0){
    sub.stage2 <- which(result == 0.5)
    data.new <- rbind(data, matrix(NA, nrow = n2, ncol = J*K))
    ## generate data only for promising subgroups
    for (i in sub.stage2){
      data.new[((n1+1):(n1+n2)),i] <- rnorm(n2, mean = mu[i], sd = sigma[i])
    }
    data <- data.new
    ## Calculate JSD for all possible subgroup divisions
    dis <- list()
    for (i in 1:length(G.futile)){
      matrixfutile <- G.futile[[i]]
      if (i == 1 | i == J*K-1){
        ## When the size of G.futile == 1, there is only one possible situation and we just need to calculate one JSD
        temp.futile <- matrixfutile
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        dis[[i]] <- calculate.Dnew(data, temp.futile, temp.eff)
      } else {
        result.dis <- c()
        for (j in 1:nrow(matrixfutile)){
          ## When the size of G.futile > 1, there is more than one possible situations and we need to calculate nrow(matrixfutile) JSDs
          temp.futile <- matrixfutile[j,]
          temp.eff <- setdiff(seq(1,J*K), temp.futile)
          result.dis[j] <- calculate.Dnew(data, temp.futile, temp.eff)
        }
        dis[[i]] <- result.dis
      }
    }
    ## Target the optimal subgroup division with JSD.max
    dis <- as.numeric(unlist(dis))
    dis.max <- which(dis == max(dis))[1] ## Obtain the serial number of the optimal subgroup division
    for (i in 1:length(G.futile)){ ## Find the high-efficacy(cluster.eff) and low-efficacy(cluster.futile) subset by recursion
      matrixfutile <- G.futile[[i]]
      if (i == 1 | i == J*K-1){
        temp.futile <- matrixfutile
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        dis.max <- dis.max - 1
        if (dis.max == 0){cluster.futile <- temp.futile; cluster.eff <- temp.eff; break}
      }
      else {
        for (j in 1:nrow(matrixfutile)){
          temp.futile <- matrixfutile[j,]
          temp.eff <- setdiff(seq(1,J*K), temp.futile)
          dis.max <- dis.max - 1
          if (dis.max == 0){cluster.futile <- temp.futile; cluster.eff <- temp.eff; break}
        }
      }
    }
    ## inference
    result <- c()
    result0 <- inference.2nd(data, cluster.futile, cluster.eff, phi.e) ## Call the inference.2nd function in "function1.R"
    result[cluster.eff] <- result0[[1]]
    result[cluster.futile] <- result0[[2]]
    ## Make the inferences conform to the monotonicity assumption
    result0 <- matrix(result, nrow = K, byrow = TRUE)
    mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
    if (sum(result) > 0){
      one <- which(result == 1)
      for (i in one){
        place <- which(mm == i, arr.ind = TRUE) ## Target the subgroup with positive result
        result0[place[1]:K, place[2]:J] <- 1
      }
      result <- as.numeric(t(result0))
    }
    N.actual <- apply(data, 2, function(x){sum(is.na(x)==FALSE)}) ## actually enrolled sample size
  }
  ## summary
  FWER <- as.numeric(sum(result[which(mu <= 0)]) > 0)  ## Family-wise type I error rate
  Power <- as.numeric(sum(result[which(mu >= 1)]) == sum(mu >= 1))  ## Conjunctive power
  if (all(mu >= 1)) {FWER <- NA}
  if (all(mu <= 0)) {Power <- NA}
  
  EN <- sum(N.actual)  ## expected sample size
  return(c(FWER, Power, EN))
}