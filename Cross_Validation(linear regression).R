##################################################
## Implementing CV methods for a simulated dataset
## We will generate X and y and generating new independent
## X.new and y.new, we will calculate the "true" prediction error
## Then we implement the LOOCV, and K-fold CV
##################################################

set.seed(10)

n <- 100
p <- 50
sigma2.star <- 4
beta.star <- rnorm(p, mean = 2)
beta.star  # to output

# repeat 500 times
B <- 5e2
foo <- 0

# Setting size of test data and training data
test.size <- floor(n/5) 

# we will use the following code to make our K-fold splits
permutation <- sample(1:n, replace = FALSE) # random permutation of 1:n
K <- 10
# Making a list of indices for each split
(test.index <- split(permutation, rep(1:K, length = n, each = n/K)))

CV.error <- matrix(0, nrow = B, ncol = 4)
time <- matrix(0, nrow = B, ncol = 3)
colnames(CV.error) <- c( "Truth","LOOCV", "10-fold", "5-fold")

for(b in 1:B)
{
  # code takes a while, hence printing for feedback
  if(b %% 100 == 0) print(b) 
  
  # Generate new d
  # Making design matrix, first column is 1
  X <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))
  
  # Generating response
  y <-  X %*% beta.star + rnorm(n, mean = 0, sd = sqrt(sigma2.star))
  
  beta.mle <- solve( t(X) %*% X) %*% t(X) %*% y
  
  # independent new X and y 
  X.new <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))
  y.new <-  X.new %*% beta.star + rnorm(n, mean = 0, sd = sqrt(sigma2.star))
  
  CV.error[b, 1] <-  mean( (y.new - X.new %*% beta.mle)^2 )
  
  
  
  ################################################
  # Leave-one-out Cross-validation
  
  time0 <- proc.time()[3]
  foo2 <- 0
  for(i in 1:n)
  {
    # Making training data
    X.train <- X[-i,] # removing ith X
    y.train <- y[-i]  #removing ith y
    
    # fitting model for training data
    beta.train <- solve(t(X.train) %*% X.train) %*% t(X.train) %*% y.train
    
    # test error
    foo2 <- foo2 + (y[i] - X[i,] %*% beta.train)^2
  }
  CV.error[b, 2] <- foo2/n
  time[b,1] <- proc.time()[3] - time0
  
  
  ################################################
  # 10-fold Cross-validation
  permutation <- sample(1:n, replace = FALSE)
  K <- 10
  
  # Making a list of indices for each split
  test.index <- split(permutation, rep(1:K, length = n, each = n/K))
  
  time0 <- proc.time()[3]
  foo3 <- 0
  for(k in 1:K)
  {
    X.train <- X[-test.index[[k]], ]
    y.train <- y[-test.index[[k]]]
    X.test <- X[test.index[[k]], ]
    y.test <- y[test.index[[k]]]
    
    beta.train <- solve(t(X.train) %*% X.train) %*% t(X.train) %*% y.train
    
    foo3 <- foo3 + sum( (y.test - X.test %*% beta.train)^2)
  }
  CV.error[b, 3] <- foo3/n
  time[b,2] <- proc.time()[3] - time0
  
  ################################################
  # 5-fold Cross-validation
  
  # Making a permutation from 1:n
  permutation <- sample(1:n, replace = FALSE)
  K <- 5
  test.index <- split(permutation, rep(1:K, length = n, each = n/K))
  
  time0 <- proc.time()[3]
  foo4 <- 0
  for(k in 1:K)
  {
    X.train <- X[-test.index[[k]], ]
    y.train <- y[-test.index[[k]]]
    X.test <- X[test.index[[k]], ]
    y.test <- y[test.index[[k]]]
    
    beta.train <- solve(t(X.train) %*% X.train) %*% t(X.train) %*% y.train
    
    foo4 <- foo4 + sum( (y.test - X.test %*% beta.train)^2)
  }
  CV.error[b, 4] <- foo4/n
  time[b,3] <- proc.time()[3] - time0
  ################################################
}

# Final estimates
colMeans(CV.error)

#LOOCV is most expensive
colMeans(time)

# LOOCV is the closest
plot(density(CV.error[,1]), xlim = range(CV.error), ylim = c(0, .30), main = "Estimated density of CV error")
lines(density(CV.error[,2]), col = "red")
lines(density(CV.error[,3]), col = "orange")
lines(density(CV.error[,4]), col = "green")
legend("topright", legend = colnames(CV.error), col = c("black", "red", "orange", "green"), lty = 1)

