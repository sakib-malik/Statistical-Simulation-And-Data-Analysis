set.seed(10)
library(mvtnorm)
library(ellipse)

data(faithful)
plot(faithful)

# calculate log-likelihood
# of mixture of multivariate normal
log_like <- function(X, pi.list, mu.list, Sigma.list, C)
{
  foo <- 0
  for(c in 1:C)
  {
    foo <- foo + pi.list[[c]]*dmvnorm(X, mean = mu.list[[c]], sigma = Sigma.list[[c]])
  }
  return(sum(log(foo)))
}

# X is the data
# C is the number of clusters
mvnEM <- function(X, C = 2, tol = 1e-3, maxit = 1e3)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  ######## Starting values ###################
  ## pi are equally split over C
  pi.list <- rep(1/C, C)
  
  mu <- list()
  Sigma <- list()
  
  # The means for each C cannot be the same, 
  # since then the three distributions overlap
  # Hence adding random noise to colMeans(X)
  for(i in 1:C)
  {
    mu[[i]] <-  rnorm(p, sd = 3) + colMeans(X)
    Sigma[[i]] <- var(X)
  }
  # Choosing good starting values is important since
  # The GMM likelihood is not concave, so the algorithm
  # may converge to a local optima.
  
  
  ######## EM algorithm steps ###################
  
  iter <- 0
  diff <- 100
  old.mu <- mu
  old.Sigma <- Sigma
  old.pi <- pi.list
  
  Ep <- matrix(0, nrow = n, ncol = C)  # gamma_{i,c}
  save.loglike <- 0
  while((diff > tol) && (iter < maxit) )
  {
    iter <- iter + 1
    ## E step: find gammas
    for(c in 1:C)
    {
      Ep[ ,c] <- pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
    }
    Ep <- Ep/rowSums(Ep)
    
    ### M-step
    pi.list <- colMeans(Ep)
    for(c in 1:C)
    {
      mu[[c]] <- colSums(Ep[ ,c] * X )/sum(Ep[,c])
    }
    
    for(c in 1:C)
    {
      foo <- 0
      for(i in 1:n)
      {
        foo <- foo + (X[i, ] - mu[[c]]) %*% t(X[i, ] - mu[[c]]) * Ep[i,c] 
      }
      Sigma[[c]] <- foo/sum(Ep[,c])
      
      # Below is to ensure the estimator is positive definite
      # otherwise next iteration gamma_i,c,k cannot be calculated
      Sigma[[c]] <- Sigma[[c]] + diag(1e-5, p)
    }
    
    save.loglike <- c(save.loglike, log_like(X = X, pi.list = pi.list, mu.list = mu, Sigma.list = Sigma,  C = C))
    # Difference in the log-likelihoods as the difference criterion
    diff <- abs(save.loglike[iter+1] - save.loglike[iter])
    
    old.mu <- mu
    old.Sigma <- Sigma
    old.pi <- pi.list
  }
  
  # Final allocation updates
  for(c in 1:C)
  {
    Ep[ ,c] <- pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
  }
  Ep <- Ep/rowSums(Ep)
  
  return(list("pi" = pi.list, "mu" = mu, "Sigma" = Sigma, "Ep" = Ep, "log.like" = tail(save.loglike,1)))
}



# We now fit this to the 
data(faithful)
C <- 2
X <- as.matrix(faithful)
class2.1 <- mvnEM(X = X, C = C)
class2.2 <- mvnEM(X = X, C = C)
class2.3 <- mvnEM(X = X, C = C)
class2.4 <- mvnEM(X = X, C = C)

# model 1,2,3 converge to a local maxima with
# log-likelihood values smaller
# model 4 converges to a point with the highest
# log-likelihood value. Hence choose model 4

par(mfrow = c(2,2))
for(i in 1:4)
{
  model <- get(paste("class2.",i, sep = ""))
  allot <- apply(model$Ep, 1, which.max)  ## Final allotment of classification
  plot(X[,1], X[,2], col = allot, pch = 16,
       main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment

  ell <- list()
  for(c in 1:C)
  {
    ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
    lines(ell[[c]], col = c)
  }  
}
# Model 4 has the largest log-likelihood, so that's the chosen model


C <- 3
class3.1 <- mvnEM(X = X, C = C)
class3.2 <- mvnEM(X = X, C = C)
class3.3 <- mvnEM(X = X, C = C)
class3.4 <- mvnEM(X = X, C = C)
class3.5 <- mvnEM(X = X, C = C)
class3.6 <- mvnEM(X = X, C = C)


par(mfrow = c(2,3))
for(i in 1:6)
{
  model <- get(paste("class3.",i, sep = ""))
  allot <- apply(model$Ep, 1, which.max)  ## Final allotment of classification
  plot(X[,1], X[,2], col = allot, pch = 16,
       main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment
  
  ell <- list()
  for(c in 1:C)
  {
    ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
    lines(ell[[c]], col = c)
  }  
}

# Model 5 has the largest log-likelihood
# so that is chosen. Here I present the final
# allocation
par(mfrow = c(1,1))
model <- class3.5
allot <- apply(model$Ep, 1, which.max)  ## Final allotment of classification
plot(X[,1], X[,2], col = allot, pch = 16,
     main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment

ell <- list()
for(c in 1:C)
{
  ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
  lines(ell[[c]], col = c)
}  

