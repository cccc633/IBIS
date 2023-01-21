Independent <- function(trial, K, J, n, mu, sigma, phi){
  ## This function is used to generate one-stage simulated trial with independent analysis ##
  ## Biomarker1 has K levels; Biomarker2 has J levels
  ## n, sample size
  ## mu, mean treatment effects for each subgroup
  ## sigma, the standard deviation of treatment effect for each subgroup
  ## phi, the statistical threshold
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## inference
  result <- c()
  for (i in 1:(J*K)){
    result[i] <- t.test(data[,i])$statistic > phi
  }
  
  ## Make the inferences conform to the monotonicity assumption
  result0 <- matrix(result, nrow = K, byrow = TRUE)
  mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
  if (sum(result) > 0){  
    one <- which(result == 1)
    for (i in one){
      place <- which(mm == i, arr.ind = TRUE)  ## Target the subgroup with positive result
      result0[place[1]:K, place[2]:J] <- 1
    }
    result <- as.numeric(t(result0))
  }
  FWER <- as.numeric(sum(result[which(mu <= 0)]) > 0)  ## Family-wise type I error rate
  ConPower <- 1 - as.numeric(sum(result[which(mu >= 1)]) < sum(mu >= 1))  ## Conjunctive power
  if (all(mu >= 1)) {FWER <- NA}
  if (all(mu <= 0)) {ConPower <- NA}
  return(list(result, FWER, ConPower))
}

BHM <- function(trial, K, J, n, mu, sigma, phi){
  ## This function is used to generate one-stage simulated trial with BHM analysis ##
  ## Biomarker1 has K levels; Biomarker2 has J levels
  ## n, sample size
  ## mu, mean treatment effects for each subgroup
  ## sigma, the standard deviation of treatment effect for each subgroup
  ## phi, the statistical threshold
  library(R2jags)
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## inference
  jags.data <- list("mean"=colMeans(data), "Ng"=J*K, "tau"=n/sd(data)^2)
  jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100)
  jags.fit <- jags.model(file = "BHM.txt", data = jags.data, inits = jags.init,
                         n.adapt=1000, n.chains=1, quiet=TRUE)
  output <- coda.samples(jags.fit, variable.names = c("mu"), n.iter=7500, thin = 3, progress.bar = "none")
  output <- as.data.frame(output[[1]])
  ## When x is the posterior samplings for theta, mean(x>0)/mean(x<=0) denotes Pr(theta>0)/Pr(theta<0)
  result <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)}) > phi
  
  ## Make the inferences conform to the monotonicity assumption
  result0 <- matrix(result, nrow = K, byrow = TRUE)
  mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
  if (sum(result) > 0){  
    one <- which(result == 1)
    for (i in one){
      place <- which(mm == i, arr.ind = TRUE)  ## Target the subgroup with positive result
      result0[place[1]:K, place[2]:J] <- 1
    }
    result <- as.numeric(t(result0))
  }
  FWER <- as.numeric(sum(result[which(mu <= 0)]) > 0)  ## Family-wise type I error rate
  ConPower <- 1 - as.numeric(sum(result[which(mu >= 1)]) < sum(mu >= 1))  ## Conjunctive power
  if (all(mu >= 1)) {FWER <- NA}
  if (all(mu <= 0)) {ConPower <- NA}
  return(list(result, FWER, ConPower))
}

IBIS <- function(trial, K, J, n, mu, sigma, G.futile, phi){
  ## This function is used to generate one-stage simulated trial with IBIS analysis ##
  ## Biomarker1 has K levels; Biomarker2 has J levels
  ## n, sample size
  ## mu, mean treatment effects for each subgroup
  ## sigma, the standard deviation of treatment effect for each subgroup
  ## G.futile, all possible low-efficacy subsets (obtained by “Divide34” function in “function1.R") 
  ## phi, the statistical threshold
  source("function1.R")
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## Calculate JSD for all possible subgroup divisions
  dis <- list()
  for (i in 1:length(G.futile)){
    matrixfutile <- G.futile[[i]]
    if (i == 1 | i == J*K-1){
      ## When the size of G.futile == 1, there is only one possible situation and we just need to calculate one JSD
      temp.futile <- matrixfutile
      temp.eff <- setdiff(seq(1,J*K), temp.futile)
      dis[[i]] <- calculate.D(n, data, temp.futile, temp.eff)
    } else {
      result.dis <- c()
      for (j in 1:nrow(matrixfutile)){
        ## When the size of G.futile > 1, there is more than one possible situations and we need to calculate nrow(matrixfutile) JSDs
        temp.futile <- matrixfutile[j,]
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        result.dis[j] <- calculate.D(n, data, temp.futile, temp.eff)
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
      result.t <- c()
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
  result0 <- inference(n, data, cluster.futile, cluster.eff, phi) ## Call the inference function in "function1.R"
  result[cluster.eff] <- result0[[1]]
  result[cluster.futile] <- result0[[2]]
  ## Make the inferences conform to the monotonicity assumption
  result0 <- matrix(result, nrow = K, byrow = TRUE)
  mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
  if (sum(result) > 0){  
    one <- which(result == 1)
    for (i in one){
      place <- which(mm == i, arr.ind = TRUE)  ## Target the subgroup with positive result
      result0[place[1]:K, place[2]:J] <- 1
    }
    result <- as.numeric(t(result0))
  }
  FWER <- as.numeric(sum(result[which(mu <= 0)]) > 0)  ## Family-wise type I error rate
  ConPower <- 1 - as.numeric(sum(result[which(mu >= 1)]) < sum(mu >= 1))  ## Conjunctive power
  if (all(mu >= 1)) {FWER <- NA}
  if (all(mu <= 0)) {ConPower <- NA}
  return(list(result, FWER, ConPower))
}

Freq <- function(trial, K, J, n, mu, sigma, G.futile, phi){
  ## This function is used to generate one-stage simulated trial with Freq analysis ##
  ## Biomarker1 has K levels; Biomarker2 has J levels
  ## n, sample size
  ## mu, mean treatment effects for each subgroup
  ## sigma, the standard deviation of treatment effect for each subgroup
  ## G.futile, all possible low-efficacy subsets (obtained by “Divide34” function in “function1.R") 
  ## phi, the statistical threshold
  source("function1.R")
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## Calculate t statistic for all possible subgroup divisions
  t.statistic <- list()
  for (i in 1:length(G.futile)){
    matrixfutile <- G.futile[[i]]
    if (i == 1 | i == J*K-1){
      ## When the size of G.futile == 1, there is only one possible situation and we just need to calculate one t statistic
      temp.futile <- matrixfutile
      temp.eff <- setdiff(seq(1,J*K), temp.futile)
      t.statistic[[i]] <- t.test(data[,temp.eff])$statistic
    } else {
      ## When the size of G.futile > 1, there is more than one possible situations and we need to calculate nrow(matrixfutile) t statistics
      result.t <- c()
      for (j in 1:nrow(matrixfutile)){
        temp.futile <- matrixfutile[j,]
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        result.t[j] <- t.test(data[,temp.eff])$statistic
      }
      t.statistic[[i]] <- result.t
    }
  }
  t.statistic[[length(G.futile)+1]] <- t.test(data)$statistic
  
  ## Target the optimal subgroup division with maximum t statistic
  t.statistic <- as.numeric(unlist(t.statistic))
  t.max <- which(t.statistic == max(t.statistic))[1] ## Obtain the serial number of the optimal subgroup division
  cluster.futile <- c(); cluster.eff <- c()
  if (t.max == length(t.statistic)){
    cluster.eff <- seq(1,J*K);
  } else {
    for (i in 1:length(G.futile)){ ## Find the high-efficacy(cluster.eff) and low-efficacy(cluster.futile) subset by recursion
      matrixfutile <- G.futile[[i]]
      if (i == 1 | i == J*K-1){
        temp.futile <- matrixfutile
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        t.max <- t.max - 1
        if (t.max == 0){cluster.futile <- temp.futile; cluster.eff <- temp.eff; break}
      }
      else {
        for (j in 1:nrow(matrixfutile)){
          temp.futile <- matrixfutile[j,]
          temp.eff <- setdiff(seq(1,J*K), temp.futile)
          t.max <- t.max - 1
          if (t.max == 0){cluster.futile <- temp.futile; cluster.eff <- temp.eff; break}
        }
      }
    }
  }
  ## inference for two subsets
  result <- rep(0, J*K)
  if (t.test(data[,cluster.eff])$statistic > phi){result[cluster.eff] <- 1}
  if (is.null(cluster.futile) == FALSE){
    if (t.test(data[,cluster.futile])$statistic > phi){result[cluster.futile] <- 1}
  }
  FWER <- as.numeric(sum(result[which(mu <= 0)]) > 0)  ## Family-wise type I error rate
  ConPower <- 1 - as.numeric(sum(result[which(mu >= 1)]) < sum(mu >= 1))  ## Conjunctive power
  if (all(mu >= 1)) {FWER <- NA}
  if (all(mu <= 0)) {ConPower <- NA}
  return(list(result, FWER, ConPower))
}

