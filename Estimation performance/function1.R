Divide34 <- function(){
  ## This function is used to enumerate all possible subgroup divisions for 3*4 scenarios (3 rows & 4 columns) ##
  
  ## Run 1000000 (a sufficient large number) simulations to get all possible subgroup orderings
  nsims <- 1000000
  y.store <- mat.or.vec(nsims,12)  ## The number of subgroups is 12
  orders.store <- mat.or.vec(nsims,12)
  ## Subgroups label
  ##  1  ##  2   ##  3   ## 4 ##
  ##  5  ##  6   ##  7   ## 8 ##
  ##  9  ##  10  ##  11  ## 12 ## 
  ##  From top to bottom, from left to right, the efficacy is increasing (marginal monotonic assumption)
  y <- rep(NA, 12)  ## y indicates the magnitude of efficacy
  y[1] <- 0  ## Subgroup 1 is the one with minimum efficacy
  y[12] <- 1 ## Subgroup 12 is the one with maximum efficacy
  set.seed(100)
  for (i in 1:nsims){
    ## Assign random values for each subgroup based on marginal monotonic assumption
    y[2] <- runif(1,min=y[1])
    y[3] <- runif(1,min=y[2])
    y[4] <- runif(1,min=y[3])
    
    y[5] <- runif(1,min=y[1])
    y[6] <- runif(1,min=max(y[5],y[2]))
    y[7] <- runif(1,min=max(y[6],y[3]))
    y[8] <- runif(1,min=max(y[7],y[4]))
    
    y[9]  <- runif(1,min=y[5])
    y[10] <- runif(1,min=max(y[9],y[6]))
    y[11] <- runif(1,min=max(y[10],y[7]))
    
    y.store[i,] <- y
  }
  ## Get all possible subgroup orderings
  for (i in 1:nsims){
    orders.store[i,] <- order(y.store[i,]) 
  }
  ## Eliminate duplicates 
  orderings34 <- unique(orders.store)
  
  ## Sort out the all possible low-efficacy subsets
  G.futile <- list()
  G.futile[[1]] <- 1 ## The low-efficacy subset include only subgroup 1
  G.futile[[11]] <- seq(1,11) ## The low-efficacy subset include all subgroups except subgroup 12
  for (i in 2:10){
    ## Find low-efficacy subset with 2 to 10 subgroups from orderings34
    noeff0 <- orderings34[,1:i]
    noeff <- t(apply(noeff0, 1, sort))
    ## Eliminate duplicates
    noeff <- unique(noeff)
    
    G.futile[[i]] <- noeff
  }
  return(G.futile)
}
calculate.D <- function(n, data, temp.futile, temp.eff){
  ## This function is used to calculate Jensen–Shannon divergence (JSD) ##
  ## n, sample size; data, trial data;
  ## temp.futile indicates which subgroups are in the low-efficacy subset
  ## temp.eff indicates which subgroups are in the high-efficacy subset
  library(R2jags)
  ## Get the posteriors of average treatment effects for both high-efficacy and low-efficacy subsets
  ## IF one subset contains only 1 subgroup, BHM cannot be applied. So there are 3 IF statements as follows.
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
  
  ## Calculate JSD with numerical calculations
  combine <- c(post.futile, post.eff)
  interval <- quantile(combine, probs = seq(0, 1, by=0.01)) ## Divide (0,1) into many intervals
  P <- c(); Q <- c() 
  for (m in 1:100){
    ## Count how many posterior samples are located in each interval to approximate the distributions
    if (m < 100){
      P[m] <- length(post.futile[(post.futile >= interval[m] & post.futile < interval[m+1])])
      Q[m] <- length(post.eff[(post.eff >= interval[m] & post.eff < interval[m+1])])
    }
    else if (m == 100){
      P[m] <- length(post.futile[(post.futile >= interval[m] & post.futile <= interval[m+1])])
      Q[m] <- length(post.eff[(post.eff >= interval[m] & post.eff <= interval[m+1])])
    }
  }
  P = P + 0.00001 ## Add a small value to avoid the probability over an interval equals 0
  Q = Q + 0.00001 
  N <- 0.5 * (P/sum(P) + Q/sum(Q))
  JS <- 0.5 * (sum(P/sum(P) * log(P/sum(P)/N, base = 2)) + sum(Q/sum(Q) * log(Q/sum(Q)/N, base = 2)))
  return(JS)
}
estimate <- function(n, data, temp.futile, temp.eff){
  ## This function is used to get the estimates with the proposed design ##
  ## n, sample size; data, trial data;
  ## temp.futile indicates which subgroups are in the low-efficacy subset
  ## temp.eff indicates which subgroups are in the high-efficacy subset
  library(R2jags)
  ## Get the posteriors of the treatment effect parameter for each subgroup
  ## IF one subset contains only 1 subgroup, BHM cannot be applied. So there are 3 IF statements as follows.
  if (length(temp.futile)==1){
    jags.data <- list("meanf"=mean(data[,temp.futile]),
                      "meane"=colMeans(data[,temp.eff]), "Ne"=length(temp.eff), 
                      "tauf"=n/sd(data[,temp.futile])^2, "taue"=n/sd(data[,temp.eff])^2)
    jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100) 
    jags.fit <- jags.model(file = "BHMcomb1f.txt", data = jags.data, inits = jags.init,
                           n.adapt=1000, n.chains=1, quiet=TRUE)
    post <- coda.samples(jags.fit, variable.names = c("mue", "muf0"), n.iter=7500, thin = 3, progress.bar = "none")
    mean.eff <- summary(post)$statistics[1:(ncol(data)-1),1]; ## mean denotes posterior mean
    lb.eff <- summary(post)$quantiles[1:(ncol(data)-1),1];    ## lb denotes the lower bound of 95% CI
    ub.eff <- summary(post)$quantiles[1:(ncol(data)-1),5];    ## ub denotes the upper bound of 95% CI
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



