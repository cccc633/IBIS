##### Defining the classification for 4*4 scenarios' efficacy #####
Divide44 <- function(){
  nsims<-1000000
  y.store<-mat.or.vec(nsims,16)
  orders.store<-mat.or.vec(nsims,16)
  y<-mat.or.vec(16,1)
  y[1]<-0
  y[16]<-1
  
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
    y[12]<-runif(1,min=max(y[11],y[8]))
    
    y[13]<-runif(1,min=y[9])
    y[14]<-runif(1,min=max(y[13],y[10]))
    y[15]<-runif(1,min=max(y[14],y[11]))
    y.store[i,]<-y
  }
  
  
  for (i in 1:nsims){
    orders.store[i,]<-order(y.store[i,]) 
  }
  
  orderings44<-unique(orders.store)
  
  count <- 0
  G.futile <- list()
  G.futile[[1]] <- 1
  for (i in 2:14){
    noeff0 <- orderings44[,1:i]
    
    noeff <- t(apply(noeff0, 1, sort))
    
    noeff <- unique(noeff)
    G.futile[[i]] <- noeff
    
    count <- count + nrow(noeff)
  }
  G.futile[[15]] <- seq(1,15)
  return(G.futile)
}
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
##### Estimation method of the proposed design #####
estimate <- function(n, data, temp.futile, temp.eff){
  library(R2jags)
  if (length(temp.futile)==1){
    jags.data <- list("meanf"=mean(data[,temp.futile]),
                      "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                      "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
    jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
    jags.fit <- jags.model(file = "BHMcomb1f.txt", data = jags.data, inits = jags.init,
                           n.adapt=1000, n.chains=1, quiet=TRUE)
    post <- coda.samples(jags.fit, variable.names = c("mue", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
    mean.eff <- summary(post)$statistics[1:(ncol(data)-1),1];
    lb.eff <- summary(post)$quantiles[1:(ncol(data)-1),1];
    ub.eff <- summary(post)$quantiles[1:(ncol(data)-1),5];
    mean.futile <- summary(post)$statistics[ncol(data),1]; 
    lb.futile <- summary(post)$quantiles[ncol(data),1];
    ub.futile <- summary(post)$quantiles[ncol(data),5];
  } 
  if (length(temp.eff)==1){
    jags.data <- list("meane"=mean(data[,temp.eff]),
                      "meanf"=colMeans(data[,temp.futile]), "Nf"=length(temp.futile), 
                      "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
    jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
    jags.fit <- jags.model(file = "BHMcomb1e.txt", data = jags.data, inits = jags.init,
                           n.adapt=1000, n.chains=1, quiet=TRUE)
    post <- coda.samples(jags.fit, variable.names = c("mue0", "muf"), n.iter=7500, thin = 3, progress.bar = "none")
    mean.eff <- summary(post)$statistics[1,1];
    lb.eff <- summary(post)$quantiles[1,1];
    ub.eff <- summary(post)$quantiles[1,5];
    mean.futile <- summary(post)$statistics[2:ncol(data),1]; 
    lb.futile <- summary(post)$quantiles[2:ncol(data),1];
    ub.futile <- summary(post)$quantiles[2:ncol(data),5];
  }
  if (length(temp.futile)>1 & length(temp.eff)>1){
    jags.data <- list("meanf"=colMeans(data[,temp.futile]), "Nf"=length(temp.futile), 
                      "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                      "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
    jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
    jags.fit <- jags.model(file = "BHMcomb.txt", data = jags.data, inits = jags.init,
                           n.adapt=1000, n.chains=1, quiet=TRUE)
    post <- coda.samples(jags.fit, variable.names = c("mue", "muf"), n.iter=7500, thin = 3, progress.bar = "none")
    mean.eff <- summary(post)$statistics[1:length(temp.eff),1];
    lb.eff <- summary(post)$quantiles[1:length(temp.eff),1];
    ub.eff <- summary(post)$quantiles[1:length(temp.eff),5];
    mean.futile <- summary(post)$statistics[(length(temp.eff)+1):ncol(data),1];
    lb.futile <- summary(post)$quantiles[(length(temp.eff)+1):ncol(data),1];
    ub.futile <- summary(post)$quantiles[(length(temp.eff)+1):ncol(data),5];
  }
  return(list(mean.futile, lb.futile, ub.futile, mean.eff, lb.eff, ub.eff))
}



