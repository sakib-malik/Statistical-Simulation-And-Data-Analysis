###########################################
## Linchpin variable sampler
## for Bayesian linear regression for cars
###########################################
BLR <- function(y, X, N = 1e4, a = 1, b = 1)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  A <- t(X)%*%X + diag(p)
  A.inv <- solve(A)
  sig2 <- numeric(length = N)
  beta <- matrix(0, nrow = N, ncol = p)
  
  rate.sig <- ( t(y) %*% (diag(1,n) - X %*% A.inv %*% t(X)) %*% y )/2 + b
  
  # sampling Inverge Gamma for sigma2
  sig2 <- 1 / rgamma(N, shape = n/2 + a, rate = rate.sig)    
  
  # Sampling beta from multivariate normal
  # mean + sqrt(covariance) %*% rnorm
  foo <- svd(A.inv)  #Singular values decomposition of A^{-1}
  Ainv.sqrt <- foo$u %*% diag(foo$d^(1/2)) %*% t(foo$v)
  
  for(i in 1:N)
  {
    beta[i,] <- A.inv %*% t(X) %*%y + Ainv.sqrt %*% rnorm(p, sd = sqrt(sig2[i]))   # Getting beta estimates
  }
  return(cbind(beta,sig2))
}

#
############################################################
#### Comparing Ridge and MLE estimates
############################################################
# We will generate similar data now, but with a Sigma correlation structure for the Xs
set.seed(1)
n <- 50
p <- 30
sigma2.star <- 1/2
beta.star <- rnorm(p)
beta.star  # to output

gen.data <- function() # Making a function to call multiple times
{
  ## setting Sigma matrix to the AR correlation matrix of .75
  foo <- matrix(.95, nrow = p-1, ncol = p-1)
  Sigma <- foo^(abs(col(foo)-row(foo)))
  Sigma
  
  # Finding sqrt(Sigma)
  decomp <- svd(Sigma)
  sqrt.Sig.svd <- decomp$u %*% diag(sqrt(decomp$d), p-1) %*% t(decomp$v)  # Symmetric matrix produced as well
  sqrt.Sig.svd %*% t(sqrt.Sig.svd) - Sigma 
  
  # Setting the rows of X_{-1} matrix to be N(0, Sigma) draws
  Xminus1 <- t(sqrt.Sig.svd %*% t(matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1))))
  
  # Making design matrix, first column is 1
  X <- cbind(1, Xminus1)
  
  # Generating response
  y <-  X %*% beta.star + rnorm(n, mean = 0, sd = sqrt(sigma2.star))
  return(list(X = X, y = y))
}

## Run once
my.data <- gen.data()
X <- my.data$X
y <- my.data$y

lambda <- 1
## Ridge solution
solve( t(X) %*%X + diag(lambda,p) ) %*% t(X) %*%y

## MLE solution
solve( t(X) %*%X ) %*% t(X) %*%y

## Now I will repeat 1000 times and store
## The average of (beta.hat - beta.star)^2

B <- 100
err.mle <- matrix(0, nrow = B, ncol = p)
err.ridge <- matrix(0, nrow = B, ncol = p)
err.pm <- matrix(0, nrow = B, ncol = p)

bias.mle <- matrix(0, nrow = B, ncol = p)
bias.ridge <- matrix(0, nrow = B, ncol = p)
bias.pm <- matrix(0, nrow = B, ncol = p)

for(b in 1:B)
{
  if(b%%(B/10) == 0) print(b/B)
  my.data <- gen.data()
  X <- my.data$X
  y <- my.data$y
  
  ## Ridge solution
  ridge <- solve( t(X) %*%X + diag(lambda,p) ) %*% t(X) %*%y
  
  ## MLE solution
  mle <- solve( t(X) %*%X ) %*% t(X) %*%y
  
  post <- BLR(y = y, X = X, N = 1e3)
  pm <- colMeans(post[,1:p])
  
  err.mle[b,] <-  (mle - beta.star)^2 
  err.ridge[b,] <- (ridge - beta.star)^2 
  err.pm[b,] <- (pm - beta.star)^2 
  
  bias.mle[b,] <- mle - beta.star
  bias.ridge[b,] <- ridge - beta.star
  bias.pm[b,] <- pm - beta.star
}

# MSE of Bayesian is similar to Ridge regression
boxplot(cbind("MLE" = colMeans(err.mle), "Ridge" = colMeans(err.ridge), "Post Mean" = colMeans(err.pm)), main = "Boxplot of Mean Squared Error") # Ridge actually does better in this regard.


# Clearly MLE returns less bias, so if bias is a criterion, MLE is better
boxplot(cbind("MLE" = colMeans(bias.mle), "Ridge" = colMeans(bias.ridge), "Post Mean" = colMeans(bias.pm)), main = "Boxplot of Bias over p variables")


