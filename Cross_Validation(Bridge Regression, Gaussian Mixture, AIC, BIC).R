##################################################
## Choosing lambda in ridge regression
## using LOOCV
## dataset is mtcars
##################################################

data(mtcars)
y <- mtcars$mpg
X <- cbind(1, as.matrix(mtcars[ ,-1]))

n <- dim(X)[1]
p <- dim(X)[2]
lam.vec <- c(10^(seq(-8, 8, by = .1)))
CV.error <- numeric(length = length(lam.vec))

for(l in 1:length(lam.vec))
{
  foo2 <- 0
  lam <- lam.vec[l]
  for(i in 1:n)
  {
    # Making training data
    X.train <- X[-i,] # removing ith X
    y.train <- y[-i]  #removing ith y
    
    # fitting model for training data
    beta.train <- solve(t(X.train) %*% X.train + lam*diag(p)) %*% t(X.train) %*% y.train
    
    # test error
    foo2 <- foo2 + (y[i] - X[i,] %*% beta.train)^2
  }
  CV.error[l] <- foo2/n
}

# lambda that yields minimum error is chosen
(chosen.lam <- lam.vec[which.min(CV.error)] )

beta.final <- solve(t(X) %*% X + chosen.lam*diag(p)) %*% t(X) %*% y
beta.final


##################################################
## Choosing lambda and alpha in bridge regression
## using 4-fold dataset is mtcars
##################################################

set.seed(12)
data(mtcars)
y <- mtcars$mpg
X <- cbind(1, as.matrix(mtcars[ ,-1]))

n <- dim(X)[1]
p <- dim(X)[2]

## MM algorithm for bridge regression
## for any alpha and lambda
bridge <- function(y, X, lambda, alpha, tol = 1e-8)
{
  current <- solve( t(X) %*%X + diag(lambda,p) ) %*% t(X) %*%y # start at ridge solution
  iter <- 0
  diff <- 100 
  while( (diff > tol) && iter < 1000)
  {
    iter <- iter + 1
    # M matrix diagonals
    ms <- as.vector(  lambda/2 *( abs(current))^(alpha - 2)  )
    
    # MM update -- using qr.solve for numerical stability
    update <- qr.solve(t(X) %*% X + diag(ms, p)) %*% t(X) %*% y
    
    diff <- sum( (current - update)^2 ) 
    current <- update
  }
  return(current)
}


lam.vec <- 1:10# c(10^(seq(-8, 8, by = .1)))
alpha <- seq(1,2, by =.05)

CV.error <- matrix(0, nrow = length(lam.vec), ncol = length(alpha))
permutation <- sample(1:n, replace = FALSE)
K <- 4

# Making a list of indices for each split
test.index <- split(permutation, rep(1:K, length = n, each = n/K))

for(l in 1:length(lam.vec))
{
  lam <- lam.vec[l]
  for(a in 1:length(alpha))
  {
    foo3 <- 0
    for(k in 1:K)
    {
      X.train <- X[-test.index[[k]], ]
      y.train <- y[-test.index[[k]]]
      X.test <- X[test.index[[k]], ]
      y.test <- y[test.index[[k]]]
      
      beta.train <- bridge(y.train, X.train, lambda = lam, alpha = alpha[a])
      
      foo3 <- foo3 + sum( (y.test - X.test %*% beta.train)^2)
    }
    CV.error[l,a] <- foo3/n
  }
}

# lambda that yields minimum error is chosen
ind <- which(CV.error == min(CV.error), arr.ind = TRUE)
(chosen.alpha <- alpha[ind[2]])
(chosen.lam <- lam.vec[ind[1]] )

# Final estimates
bridge(y, X, lambda = chosen.lam, alpha = chosen.alpha)







###########################################
## EM algorithm to fit mixture of Gaussians
## to multivariate data (old faithful)
##
## Cross-validation and AIC model selection
###########################################
library(mvtnorm)
library(ellipse)
set.seed(12)


# calculate negative log-likelihood
# of mixture of multivariate normal
log_like <- function(X, pi.list, mu.list, Sigma.list, C)
{
  foo <- 0
  for(c in 1:C)
  {
    foo <- foo + pi.list[[c]]*dmvnorm(X, mean = mu.list[[c]], sigma = Sigma.list[[c]])
  }
  return(-sum(log(foo)))
}

# X is the data
# C is the number of clusters
GLMMforC <- function(X, C, tol = 1e-3, maxit = 1e3)
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



# My data set
data(faithful)
X <- as.matrix(faithful)
n <- dim(X)[1]

################################################
# 5-fold Cross-validation

permutation <- sample(1:n, replace = FALSE)
K <- 5
# Uneven folds, but that is ok
test.index <- split(permutation, rep(1:K, length = n, each = n/K))

# Testing whether 2-4 classes are needed
potC <- 2:4
CV.errorLike <- numeric(length = length(potC))

# will run EM multiple times for each training data
# since EM convergence to local minima. Setting these
# reps = 7
reps <- 5
model.save <- list()
for(c in 1:length(potC))
{
  foo3 <- 0
  for(k in 1:K)
  {
    print(c(c,k))
    X.train <- X[-test.index[[k]], ]
    X.test <- X[test.index[[k]], ]
    
    for(r in 1:reps)
    {
      model.save[[r]] <- GLMMforC(X = X.train, C = potC[c])
    }
    # which ever run is the lowest negative log-like
    chosen.run <- which.min(sapply(model.save, function(t) t$log.like))
    model <- model.save[[chosen.run]]
    foo3 <- foo3 + log_like(X = X.test, pi.list = model$pi, mu.list =model$mu, Sigma.list = model$Sigma,  C = potC[c])
  }
  CV.errorLike[c] <- foo3/n
}

CV.errorLike #Lowest value is for C = 3, 

# Thus choose C = 3 classes.

######################################
##### Model selection via AIC #######
######################################

aic <- function(X, pi.list, mu.list, Sigma.list, C)
{
  nlike <- log_like(X, pi.list, mu.list, Sigma.list, C)
  rtn <- 2*nlike + 2* (6*C - 1) # No. of params  = 3*C -  1
  return(rtn)
}

bic <- function(X, pi.list, mu.list, Sigma.list, C)
{
  n <- dim(X)[1]
  nlike <- log_like(X, pi.list, mu.list, Sigma.list, C)
  rtn <- 2*nlike + log(n)* (6*C - 1) # No. of params  = 3*C -  1
  return(rtn)
}

aicLike <- numeric(length = length(potC))
bicLike <- numeric(length = length(potC))
reps <- 5
model.save <- list()
model <- list()
for(c in 1:length(potC))
{
  print(c)
  for(r in 1:reps)
  {
    model.save[[r]] <- GLMMforC(X = X, C = potC[c])
  } 
  
  chosen.run <- which.min(sapply(model.save, function(t) t$log.like))
  model[[c]] <- model.save[[chosen.run]]
  aicLike[c] <- aic(X = X, pi.list = model[[c]]$pi, mu.list =model[[c]]$mu, Sigma.list = model[[c]]$Sigma, C = potC[c])
  bicLike[c] <- bic(X = X, pi.list = model[[c]]$pi, mu.list =model[[c]]$mu, Sigma.list = model[[c]]$Sigma, C = potC[c])
}

aicLike  # Lowest is C = 4 is the best!
bicLike  # Lowest is C = 2

par(mfrow = c(1,2))
chosen <- which.min(aicLike)
allot <- apply(model[[chosen]]$Ep, 1, which.max)  ## Final allotment of classification
plot(X[,1], X[,2], col = allot, pch = 16, main = "AIC: C = 4") # plot allotment

ell <- list()
for(c in 1:potC[[chosen]])
{
  ell[[c]] <- ellipse(model[[chosen]]$Sigma[[c]], centre = as.numeric(model[[chosen]]$mu[[c]]))
  lines(ell[[c]], col = c)
}  

chosen <- which.min(bicLike)
allot <- apply(model[[chosen]]$Ep, 1, which.max)  ## Final allotment of classification
plot(X[,1], X[,2], col = allot, pch = 16, main = "BIC: C = 2") # plot allotment

ell <- list()
for(c in 1:potC[[chosen]])
{
  ell[[c]] <- ellipse(model[[chosen]]$Sigma[[c]], centre = as.numeric(model[[chosen]]$mu[[c]]))
  lines(ell[[c]], col = c)
}  


