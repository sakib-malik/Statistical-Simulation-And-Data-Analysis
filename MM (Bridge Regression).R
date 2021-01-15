################################################
## Bridge regression and the MM algorithm
## Compare for different values of alpha
################################################
set.seed(1)
n <- 100
p <- 5
beta.star <- c(0,0,0,rnorm(p-3, sd = 1))  # larger variance than exercise 4.
beta.star  # to output

# Making design matrix, first column is 1
X <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))

# Generating response
y <-  X %*% beta.star + rnorm(n, mean = 0, sd = 1)

###############
# First MLE
###############
# MLE of beta
beta.mle <- solve( t(X) %*%X ) %*% t(X) %*%y

# Bridge solutions
alpha.vec <- c(1, 1.5, 1.8)

lambda <- 10
beta.bridge <- matrix(0, nrow = p, ncol = 3)
tol <- 1e-5
for(i in 1:length(alpha.vec) ) # loop for each alpha
{
  current <- solve( t(X) %*%X + diag(lambda,p) ) %*% t(X) %*%y # start at ridge solution
  iter <- 0
  diff <- 100 
  while( (diff > tol) && iter < 1000)
  {
    iter <- iter + 1
    
    # M matrix diagonals
    ms <- as.vector(  lambda/2 *( abs(current))^(alpha.vec[i] - 2)  )
    
    # MM update -- using qr.solve for numerical stability
    update <- qr.solve(t(X) %*% X + diag(ms, p)) %*% t(X) %*% y

    diff <- sum( (current - update)^2 ) 
    current <- update
  }
  beta.bridge[,i] <- current
}

## Comparing MLE and Bridge for different alpha values
beta.ridge <- solve( t(X) %*%X + diag(lambda,p) ) %*% t(X) %*%y

cbind(beta.mle, beta.bridge, beta.ridge)
beta.star
# Ridge is event closer to 0
