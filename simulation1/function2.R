Independent <- function(trial, K, J, n, mu, sigma){
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## estimate
  mean <- c()
  CI.lb <- c()
  CI.ub <- c()
  for (i in 1:(J*K)){
    temp <- t.test(data[,i])
    mean[i] <- temp$estimate
    CI.lb[i] <- temp$conf.int[1]
    CI.ub[i] <- temp$conf.int[2]
  }
  ## calculate metrics
  Bias <- mean - mu
  MSE <- (mean - mu)^2
  Width <- CI.ub - CI.lb
  return(list(Bias, MSE, Width))
}

BHM <- function(trial, K, J, n, mu, sigma){
  library(R2jags)
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## estimate
  jags.data <- list("mean"=colMeans(data), "Ng"=J*K, "tau"=n/sd(data)^2)
  jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100)
  # jags.init <- list("mu0"=0, "tau0"=1) 
  jags.fit <- jags.model(file = "BHM.txt", data = jags.data, inits = jags.init,
                         n.adapt=1000, n.chains=1, quiet=TRUE)
  output <- coda.samples(jags.fit, variable.names = c("mu"), n.iter=7500, thin = 3, progress.bar = "none")
  temp <- summary(output)
  mean <- temp$statistics[,1]
  CI.lb <- temp$quantiles[,1]
  CI.ub <- temp$quantiles[,5]
  ## calculate metrics
  Bias <- mean - mu
  MSE <- (mean - mu)^2
  Width <- CI.ub - CI.lb
  return(list(Bias, MSE, Width))
}

IBIS <- function(trial, K, J, n, mu, sigma, G.futile){
  source("function.R")
  set.seed(trial)
  data <- matrix(NA, nrow = n, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n, mean = mu[i], sd = sigma[i])
  }
  ## Calculate JS divergence
  distance <- list()
  for (i in 1:length(G.futile)){
    matrixfutile <- G.futile[[i]]
    if (i == 1 | i == J*K-1){
      temp.futile <- matrixfutile
      temp.eff <- setdiff(seq(1,J*K), temp.futile)
      distance[[i]] <- calculate.D(n, data, temp.futile, temp.eff)
    } else {
      result.dis <- c()
      for (j in 1:nrow(matrixfutile)){
        temp.futile <- matrixfutile[j,]
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        result.dis[j] <- calculate.D(n, data, temp.futile, temp.eff)
      }
      distance[[i]] <- result.dis
    }
  }
  ## Target JS.max
  distance <- as.numeric(unlist(distance))
  dis.max <- which(distance == max(distance))[1]
  for (i in 1:length(G.futile)){
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
  temp <- estimate(n, data, cluster.futile, cluster.eff)
  mean = CI.lb = CI.ub = c()
  mean[cluster.futile] <- temp[[1]]; mean[cluster.eff] <- temp[[4]];
  CI.lb[cluster.futile] <- temp[[2]]; CI.lb[cluster.eff] <- temp[[5]];
  CI.ub[cluster.futile] <- temp[[3]]; CI.ub[cluster.eff] <- temp[[6]];
  ## calculate metrics
  Bias <- mean - mu
  MSE <- (mean - mu)^2
  Width <- CI.ub - CI.lb
  return(list(Bias, MSE, Width))
}