

################################################
## MLE for location Cauchy distribution using
## Gradient Ascent method
## We will plot the likelihood as well
################################################
set.seed(1)
mu.star <- 5  # Setting true mu
n <- 4  # sample size
X <- rt(n, df = 1) + mu.star

## Function calculates  the log-likelihood
log.like <- function(mu, X)
{
  n <- length(X)
  rtn <- -n*log(pi) - sum( log(1 + (X - mu)^2) )
  return(rtn)
}

mu.x <- seq(-10, 40, length = 1e3)  # A sequence of mu's 
ll.est <- sapply(mu.x, log.like, X)  # evaluating log-likelihood at the mus
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))  # plotting log-likelihood. Not concave, so we need to choose good starting values.


## Starting Gradient-Ascent method
tol <- 1e-5  # tolerance level for when to stop algorithm


## Returns derivate of log-likelihood
f.prime <- function(X, mu)
{
  rtn <- sum(2* (X - mu)/(1 + (X-mu)^2))  #f.prime
  return(rtn)
}

# Returns double derivative of log-likelihood.
f.double.prime <- function(X, mu)
{
  rtn <- sum( 2 * ( 2*(X-mu)^2/ (1 + (X - mu)^2)^2  - (1 + (X-mu)^2)^(-1) )  )
  return(rtn)
}

## Loop below stops when |mu_(k+1) - mu_(k)| < tol
t <- .3  # change this to 1 and see what happens to the "bad" starting value



## Loop below stops when |mu_(k+1) - mu_(k)| < tol
current <- 7  # Bad starting value
diff <- 100
iter <- 0
mu.k <- current
while( (diff > tol) && iter < 100)
{
  iter <- iter + 1
  update <- current + t*f.prime(X, current) 
  mu.k <- c(mu.k, update)
  diff <- abs(current - update)
  current <- update
}
current  # final approximation to MLE
evals <- sapply(mu.k, log.like, X)
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))
points(mu.k, evals, pch = 16, col = rgb(1,0,0, alpha = .5))


current <- median(X)  # Good starting value
diff <- 100  # inital large value for difference
iter <- 0    # counting the number of iterations

mu.k <- current
while( (diff > tol) && iter < 100)
{
  iter <- iter + 1
  update <- current + t*f.prime(X, current)  # GD update
  mu.k <- c(mu.k, update)
  diff <- abs(current - update)
  current <- update
}
current  # final approximation to MLE
evals <- sapply(mu.k, log.like, X)
points(mu.k, evals, pch = 16, col = rgb(0,0,1, alpha = .5))


current <- 19  # Worst starting value
diff <- 100
iter <- 0
mu.k <- current
while( (diff > tol) && iter < 100)
{
  iter <- iter + 1
  update <- current + t*f.prime(X, current) 
  mu.k <- c(mu.k, update)
  diff <- abs(current - update)
  current <- update
}
current  # final approximation to MLE
evals <- sapply(mu.k, log.like, X)
points(mu.k, evals, pch = 16, col = rgb(.2,.7,.1, alpha = .8))

legend("topright", legend = c("Good starting", "Bad starting", "Horrible starting"), pch = 16, col = c("blue", "red", rgb(.2,.7,.1)))




################################################
## MLE for logistic regression
## Using gradient ascent
################################################
library(mcmc) #to load a dataset
data(logit)
head(logit)  # y is response and 4 covariates
y <- logit$y
X <- as.matrix(logit[, 2:5])
p <- dim(X)[2]

f.gradient <- function(y, X, beta)
{
  beta <- matrix(beta, ncol = 1)
  pi <- exp(X %*% beta) / (1 + exp(X%*%beta))  
  rtn <- colSums(X* as.numeric(y - pi))
  return(rtn)
}


store.beta <- matrix(0, nrow = 1, ncol = p)
beta_k <- rep(0, p) # start at all 0s
t <- .1
tol <- 1e-5
iter <- 0
diff <- 100
while((diff > tol) && iter < 100)  #not too many iterations
{
  iter <- iter+1
  old <- beta_k
  beta_k = old + t* f.gradient(y = y, X= X, beta = old)
  store.beta <- rbind(store.beta, beta_k)
  diff <- sum( (beta_k - old)^2)
}
iter # number of iterations
beta_k # last estimate





