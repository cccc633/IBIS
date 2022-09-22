##### Defining the classification for 3*4 scenarios' efficacy (3 rows & 4 columns)#####
Divide34 <- function(){
  nsims<-1000000
  y.store<-mat.or.vec(nsims,12)
  orders.store<-mat.or.vec(nsims,12)
  y<-mat.or.vec(12,1)
  y[1]<-0
  y[12]<-1
  
  set.seed(100)
  for (i in 1:nsims){
    y[2]<-runif(1,min=y[1])
    y[3]<-runif(1,min=y[2])
    y[4]<-runif(1,min=y[3])
    
    y[5]<-runif(1,min=y[1])
    y[6]<-runif(1,min=max(y[5],y[2]))
    y[7]<-runif(1,min=max(y[6],y[3]))
    y[8]<-runif(1,min=max(y[7],y[4]))
    
    y[9] <-runif(1,min=y[5])
    y[10]<-runif(1,min=max(y[9],y[6]))
    y[11]<-runif(1,min=max(y[10],y[7]))
    
    y.store[i,]<-y
  }
  
  for (i in 1:nsims){
    orders.store[i,]<-order(y.store[i,]) 
  }
  
  orderings34<-unique(orders.store)
  
  G.futile <- list()
  G.futile[[1]] <- 1
  G.futile[[11]] <- seq(1,11)
  count <- 2
  for (i in 2:10){
    noeff0 <- orderings34[,1:i]
    
    noeff <- t(apply(noeff0, 1, sort))
    
    noeff <- unique(noeff)
    G.futile[[i]] <- noeff
    
    count <- count + nrow(noeff)
  }
  return(G.futile)
}

IBIS <- function(trial, K, J, n1, n2, mu, sigma, G.futile, phi.e, phi.p){
  ##### Calculate JSD #####
  calculate.D <- function(n, data, temp.futile, temp.eff){
    library(R2jags)
    if (length(temp.futile)==1){
      jags.data <- list("meanf"=mean(data[,temp.futile]),
                        "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                        "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1f.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      post <- coda.samples(jags.fit, variable.names = c("mue0", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      post.eff <- as.data.frame(post[[1]])$mue0
      post.futile <- as.data.frame(post[[1]])$muf0
    } 
    if (length(temp.eff)==1){
      jags.data <- list("meanf"=colMeans(data[,temp.futile]), "Nf"=length(temp.futile),
                        "meane"=mean(data[,temp.eff]),
                        "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1e.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      post <- coda.samples(jags.fit, variable.names = c("mue0", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      post.eff <- as.data.frame(post[[1]])$mue0
      post.futile <- as.data.frame(post[[1]])$muf0
    }
    if (length(temp.futile)>1 & length(temp.eff)>1){
      jags.data <- list("meanf"=colMeans(data[,temp.futile]), "Nf"=length(temp.futile), 
                        "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                        "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      post <- coda.samples(jags.fit, variable.names = c("mue0", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      post.eff <- as.data.frame(post[[1]])$mue0
      post.futile <- as.data.frame(post[[1]])$muf0
    }
    ## calculate JS divergence
    combine <- c(post.futile, post.eff)
    interval <- quantile(combine, probs = seq(0, 1, by=0.01))
    P <- c(); Q <- c()
    for (m in 1:100){
      if (m < 100){
        P[m] <- length(post.futile[(post.futile >= interval[m] & post.futile < interval[m+1])])
        Q[m] <- length(post.eff[(post.eff >= interval[m] & post.eff < interval[m+1])])
      }
      else if (m == 100){
        P[m] <- length(post.futile[(post.futile >= interval[m] & post.futile <= interval[m+1])])
        Q[m] <- length(post.eff[(post.eff >= interval[m] & post.eff <= interval[m+1])])
      }
    }
    P = P + 0.00001
    Q = Q + 0.00001
    N <- 0.5 * (P/sum(P) + Q/sum(Q))
    JS <- 0.5 * (sum(P/sum(P) * log(P/sum(P)/N, base = 2)) + sum(Q/sum(Q) * log(Q/sum(Q)/N, base = 2)))
    return(JS)
  }
  ##### Inference for the proposed design #####
  inference <- function(n, data, temp.futile, temp.eff, phi.e, phi.p){
    library(R2jags)
    result <- c()
    if (length(temp.futile)==1){
      jags.data <- list("meanf"=mean(data[,temp.futile]),
                        "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                        "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1f.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      output <- coda.samples(jags.fit, variable.names = c("mue", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      output <- as.data.frame(output[[1]])
      result <- rep(0, J*K)
      BF <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)})
      result[BF > phi.p] <- 0.5 ## promising
      result[BF > phi.e] <- 1   ## effective
      result.eff <- result[1:(ncol(data)-1)]
      result.futile <- result[ncol(data)]
    } 
    if (length(temp.eff)==1){
      jags.data <- list("meane"=mean(data[,temp.eff]),
                        "meanf"=colMeans(data[,temp.futile]), "Nf"=length(temp.futile), 
                        "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1e.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      output <- coda.samples(jags.fit, variable.names = c("mue0", "muf"), n.iter=7500, thin = 3, progress.bar = "none")
      output <- as.data.frame(output[[1]])
      result <- rep(0, J*K)
      BF <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)})
      result[BF > phi.p] <- 0.5 ## promising
      result[BF > phi.e] <- 1   ## effective
      result.eff <- result[1]
      result.futile <- result[2:ncol(data)]
    }
    if (length(temp.futile)>1 & length(temp.eff)>1){
      jags.data <- list("meanf"=colMeans(data[,temp.futile]), "Nf"=length(temp.futile), 
                        "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                        "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      output <- coda.samples(jags.fit, variable.names = c("mue", "muf"), n.iter=7500, thin = 3, progress.bar = "none")
      output <- as.data.frame(output[[1]])
      result <- rep(0, J*K)
      BF <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)})
      result[BF > phi.p] <- 0.5 ## promising
      result[BF > phi.e] <- 1   ## effective
      result.eff <- result[1:length(temp.eff)]
      result.futile <- result[(length(temp.eff)+1):ncol(data)]
    }
    return(list(result.eff, result.futile))
  }
  
  ##### Calculate JSD #####
  calculate.Dnew <- function(data, temp.futile, temp.eff){
    library(R2jags)
    n <- apply(data, 2, function(x){sum(is.na(x)==FALSE)})
    if (length(temp.futile)==1){
      jags.data <- list("meanf"=mean(data[,temp.futile], na.rm = TRUE),
                        "meane"=colMeans(data[,temp.eff], na.rm = TRUE), "Ne"=length(temp.eff), 
                        "tauf"=n[temp.futile]/sd(data[,temp.futile], na.rm = TRUE)^2, 
                        "taue"=n[temp.eff]/sd(data[,temp.eff], na.rm = TRUE)^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1f-f4.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      post <- coda.samples(jags.fit, variable.names = c("mue0", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      post.eff <- as.data.frame(post[[1]])$mue0
      post.futile <- as.data.frame(post[[1]])$muf0
    } 
    if (length(temp.eff)==1){
      jags.data <- list("meanf"=colMeans(data[,temp.futile], na.rm = TRUE), "Nf"=length(temp.futile),
                        "meane"=mean(data[,temp.eff], na.rm = TRUE),
                        "tauf"=n[temp.futile]/sd(data[,temp.futile], na.rm = TRUE)^2, 
                        "taue"=n[temp.eff]/sd(data[,temp.eff], na.rm = TRUE)^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1e-f4.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      post <- coda.samples(jags.fit, variable.names = c("mue0", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      post.eff <- as.data.frame(post[[1]])$mue0
      post.futile <- as.data.frame(post[[1]])$muf0
    }
    if (length(temp.futile)>1 & length(temp.eff)>1){
      jags.data <- list("meanf"=colMeans(data[,temp.futile], na.rm = TRUE), "Nf"=length(temp.futile), 
                        "meane"=colMeans(data[,temp.eff], na.rm = TRUE), "Ne"=length(temp.eff), 
                        "tauf"=n[temp.futile]/sd(data[,temp.futile], na.rm = TRUE)^2, 
                        "taue"=n[temp.eff]/sd(data[,temp.eff], na.rm = TRUE)^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb-f4.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      post <- coda.samples(jags.fit, variable.names = c("mue0", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      post.eff <- as.data.frame(post[[1]])$mue0
      post.futile <- as.data.frame(post[[1]])$muf0
    }
    ## calculate JS divergence
    combine <- c(post.futile, post.eff)
    interval <- quantile(combine, probs = seq(0, 1, by=0.01))
    P <- c(); Q <- c()
    for (m in 1:100){
      if (m < 100){
        P[m] <- length(post.futile[(post.futile >= interval[m] & post.futile < interval[m+1])])
        Q[m] <- length(post.eff[(post.eff >= interval[m] & post.eff < interval[m+1])])
      }
      else if (m == 100){
        P[m] <- length(post.futile[(post.futile >= interval[m] & post.futile <= interval[m+1])])
        Q[m] <- length(post.eff[(post.eff >= interval[m] & post.eff <= interval[m+1])])
      }
    }
    P = P + 0.00001
    Q = Q + 0.00001
    N <- 0.5 * (P/sum(P) + Q/sum(Q))
    JS <- 0.5 * (sum(P/sum(P) * log(P/sum(P)/N, base = 2)) + sum(Q/sum(Q) * log(Q/sum(Q)/N, base = 2)))
    return(JS)
  }
  ##### Inference for the proposed design #####
  inference.new <- function(data, temp.futile, temp.eff, phi.e){
    library(R2jags)
    n <- apply(data, 2, function(x){sum(is.na(x)==FALSE)})
    result <- c()
    if (length(temp.futile)==1){
      jags.data <- list("meanf"=mean(data[,temp.futile], na.rm = TRUE),
                        "meane"=colMeans(data[,temp.eff], na.rm = TRUE), "Ne"=length(temp.eff), 
                        "tauf"=n[temp.futile]/sd(data[,temp.futile], na.rm = TRUE)^2, 
                        "taue"=n[temp.eff]/sd(data[,temp.eff], na.rm = TRUE)^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1f-f4.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      output <- coda.samples(jags.fit, variable.names = c("mue", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
      output <- as.data.frame(output[[1]])
      result <- rep(0, J*K)
      BF <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)})
      result[BF > phi.e] <- 1   ## effective
      result.eff <- result[1:(ncol(data)-1)]
      result.futile <- result[ncol(data)]
    } 
    if (length(temp.eff)==1){
      jags.data <- list("meane"=mean(data[,temp.eff], na.rm = TRUE),
                        "meanf"=colMeans(data[,temp.futile], na.rm = TRUE), "Nf"=length(temp.futile), 
                        "tauf"=n[temp.futile]/sd(data[,temp.futile], na.rm = TRUE)^2, 
                        "taue"=n[temp.eff]/sd(data[,temp.eff], na.rm = TRUE)^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb1e-f4.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      output <- coda.samples(jags.fit, variable.names = c("mue0", "muf"), n.iter=7500, thin = 3, progress.bar = "none")
      output <- as.data.frame(output[[1]])
      result <- rep(0, J*K)
      BF <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)})
      result[BF > phi.e] <- 1   ## effective
      result.eff <- result[1]
      result.futile <- result[2:ncol(data)]
    }
    if (length(temp.futile)>1 & length(temp.eff)>1){
      jags.data <- list("meanf"=colMeans(data[,temp.futile], na.rm = TRUE), "Nf"=length(temp.futile), 
                        "meane"=colMeans(data[,temp.eff], na.rm = TRUE), "Ne"=length(temp.eff), 
                        "tauf"=n[temp.futile]/sd(data[,temp.futile], na.rm = TRUE)^2, 
                        "taue"=n[temp.eff]/sd(data[,temp.eff], na.rm = TRUE)^2)
      jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
      jags.fit <- jags.model(file = "BHMcomb-f4.txt", data = jags.data, inits = jags.init,
                             n.adapt=1000, n.chains=1, quiet=TRUE)
      output <- coda.samples(jags.fit, variable.names = c("mue", "muf"), n.iter=7500, thin = 3, progress.bar = "none")
      output <- as.data.frame(output[[1]])
      result <- rep(0, J*K)
      BF <- apply(output, 2, function(x){mean(x>0)/mean(x<=0)})
      result[BF > phi.e] <- 1   ## effective
      result.eff <- result[1:length(temp.eff)]
      result.futile <- result[(length(temp.eff)+1):ncol(data)]
    }
    return(list(result.eff, result.futile))
  }
  
  ##### stage 1 #####
  set.seed(trial)
  data <- matrix(NA, nrow = n1, ncol = J*K)
  ## generate data
  for (i in 1:(J*K)){
    data[,i] <- rnorm(n1, mean = mu[i], sd = sigma[i])
  }
  ## Calculate JS divergence
  dis <- list()
  for (i in 1:length(G.futile)){
    matrixfutile <- G.futile[[i]]
    if (i == 1 | i == J*K-1){
      temp.futile <- matrixfutile
      temp.eff <- setdiff(seq(1,J*K), temp.futile)
      dis[[i]] <- calculate.D(n1, data, temp.futile, temp.eff)
    } else {
      result.dis <- c()
      for (j in 1:nrow(matrixfutile)){
        temp.futile <- matrixfutile[j,]
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        result.dis[j] <- calculate.D(n1, data, temp.futile, temp.eff)
      }
      dis[[i]] <- result.dis
    }
  }
  ## Target JS.max
  dis <- as.numeric(unlist(dis))
  dis.max <- which(dis == max(dis))[1]
  for (i in 1:length(G.futile)){
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
  result0 <- inference(n1, data, cluster.futile, cluster.eff, phi.e, phi.p)
  result[cluster.eff] <- result0[[1]]
  result[cluster.futile] <- result0[[2]]
  ## fit monotonicity hypothesis
  result0 <- matrix(result, nrow = K, byrow = TRUE)
  mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
  if (sum(result) > 0){  
    half <- which(result == 0.5)
    for (i in half){
      place <- which(mm == i, arr.ind = TRUE)
      result0[place[1]:K, place[2]:J] <- 0.5
    }
    one <- which(result == 1)
    for (i in one){
      place <- which(mm == i, arr.ind = TRUE)
      result0[place[1]:K, place[2]:J] <- 1
    }
    result <- as.numeric(t(result0))
  }
  N.actual <- rep(n1, J*K)
  ##### stage 2 #####
  if (sum(result == 0.5) > 0){
    sub.stage2 <- which(result == 0.5)
    data.new <- rbind(data, matrix(NA, nrow = n2, ncol = J*K))
    ## generate data
    for (i in sub.stage2){
      data.new[((n1+1):(n1+n2)),i] <- rnorm(n2, mean = mu[i], sd = sigma[i])
    }
    data <- data.new
    ## Calculate JS divergence
    dis <- list()
    for (i in 1:length(G.futile)){
      matrixfutile <- G.futile[[i]]
      if (i == 1 | i == J*K-1){
        temp.futile <- matrixfutile
        temp.eff <- setdiff(seq(1,J*K), temp.futile)
        dis[[i]] <- calculate.Dnew(data, temp.futile, temp.eff)
      } else {
        result.dis <- c()
        for (j in 1:nrow(matrixfutile)){
          temp.futile <- matrixfutile[j,]
          temp.eff <- setdiff(seq(1,J*K), temp.futile)
          result.dis[j] <- calculate.Dnew(data, temp.futile, temp.eff)
        }
        dis[[i]] <- result.dis
      }
    }
    ## Target JS.max
    dis <- as.numeric(unlist(dis))
    dis.max <- which(dis == max(dis))[1]
    for (i in 1:length(G.futile)){
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
    result0 <- inference.new(data, cluster.futile, cluster.eff, phi.e)
    result[cluster.eff] <- result0[[1]]
    result[cluster.futile] <- result0[[2]]
    ## fit monotonicity hypothesis
    result0 <- matrix(result, nrow = K, byrow = TRUE)
    mm <- matrix(1:(J*K), nrow = K, byrow = TRUE)
    if (sum(result) > 0){
      one <- which(result == 1)
      for (i in one){
        place <- which(mm == i, arr.ind = TRUE)
        result0[place[1]:K, place[2]:J] <- 1
      }
      result <- as.numeric(t(result0))
    }
    N.actual <- apply(data, 2, function(x){sum(is.na(x)==FALSE)})
  }
  ## summary
  FWER1 <- as.numeric(sum(result[which(mu <= 0)]) > 0)
  Power <- as.numeric(sum(result[which(mu >= 1)]) == sum(mu >= 1))
  EN <- sum(N.actual)
  return(c(FWER1, Power, EN))
}



