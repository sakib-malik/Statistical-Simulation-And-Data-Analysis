################################################
## MLE for logistic regression
## Using stochastic gradient ascent
################################################

f.gradient <- function(y, X, beta)
{
  n <- dim(X)[1]
  beta <- matrix(beta, ncol = 1)
  pi <- exp(X %*% beta) / (1 + exp(X%*%beta))  
  rtn <- colSums(X* as.numeric(y - pi))
  return(n*rtn)
}


#################################################
## The following is a general function that
## implements the regular gradient ascent
## the stochastic gradient ascent and 
## mini-batch stochastic gradient ascent
#################################################
SGA <- function(y, X, batch.size = dim(X)[1], t = .1, max.iter = dim(X)[1], adapt = FALSE)
{  
  p <- dim(X)[2]
  n <- dim(X)[1]

  # create the mini-batches
  permutation <- sample(1:n, replace = FALSE)
  K <- floor(n/batch.size)
  batch.index <- split(permutation, rep(1:K, length = n, each = n/K))

  # index for choosing the mini-batch
  count <- 1

  beta_k <- rep(0, p) # start at all 0s
  track.gradient <- matrix(0, nrow = max.iter, ncol = p)
  track.gradient[1,] <- f.gradient(y = y, X= X, beta = beta_k)/n
  
  # saving the running mean of the estimates of theta^*
  mean_beta <- rep(0,p)
  
  # tk: in case we want t_k
  tk <- t
  
  # ideally, we will have a while loop here, but
  # I have written this to always complete some max.iter steps
  for(iter in 1:max.iter)  
  {
    count <- count+1

    if(adapt) tk <- t/(sqrt(iter))  # in case t_k
    if(count %% K == 0) count <- count%%K  +1  # when all batches finish, restart the batches
    if(iter %% (max.iter/10) == 0) print(iter) #feedback

    # batch of data
    y.batch <- y[batch.index[[count]] ]
    X.batch <- matrix(X[batch.index[[count]], ], nrow = batch.size)
    
    # SGA step
    beta_k = beta_k + tk* f.gradient(y = y.batch, X = X.batch, beta = beta_k)/batch.size

    # saving overall estimates and running gradients for demonstration
    mean_beta <- (beta_k + mean_beta*(iter - 1))/(iter)
    if(batch.size == n)
    {
      est <- beta_k
    }else{
      est <-  mean_beta
    }    
    track.gradient[iter,] <- f.gradient(y = y, X = X, beta = est)/n
  }
  rtn <- list("iter" = iter, "est" = est, "grad" = track.gradient[1:iter,])
  return(rtn)
}


# Generating data for demonstration
set.seed(10)
p <- 5
n <- 1e4
X <- matrix(rnorm(n*(p-1)), nrow = n, ncol = p-1)
X <- cbind(1, X)
beta <- matrix(rnorm(p, 0, sd = 1), ncol = 1)
p <- exp(X %*% beta)/(1 + exp(X%*%beta))
y <- rbinom(n, size = 1, prob = p)

#True MLE
coeff <- glm(y ~ X -1,family = "binomial")$coeff

# Tuned value of t
ga <- SGA(y, X, batch.size = 1e4, t = .0015, max.iter = 1e3)
b1 <- SGA(y, X, batch.size = 1, t = .1, max.iter = ga$iter) 
b10 <- SGA(y, X, batch.size = 10, t = .1, max.iter = ga$iter)
b100 <- SGA(y, X, batch.size = 100, t = .1, max.iter = ga$iter) 

index <- 1:500
plot(apply(ga$grad[index,], 1, function(t) sum(abs(t))), type = 'l', ylim = c(0,max(apply(b1$grad[,], 1, function(t) sum(abs(t))))), ylab = "Complete gradient")
lines(apply(b1$grad[index,], 1, function(t) sum(abs(t))), col = "red")
lines(apply(b10$grad[index,], 1, function(t) sum(abs(t))), col = "blue")
lines(apply(b100$grad[index,], 1, function(t) sum(abs(t))), col = "orange")
legend("topright", col = c("black", "red", "blue", "orange"), legend = c("GA", "SGA", "MB-SGA-10", "MB-SGA-100"), lty = 1)

# all bad values of t
ga <- SGA(y, X, batch.size = n, t = .005, max.iter = 1e3)
b1 <- SGA(y, X, batch.size = 1, t = 1, max.iter = ga$iter) 
b10 <- SGA(y, X, batch.size = 10, t = 1, max.iter = ga$iter)
b100 <- SGA(y, X, batch.size = 100, t = 1, max.iter = ga$iter) 

index <- 1:1000
plot(apply(ga$grad[index,], 1, function(t) sum(abs(t))), type = 'l', ylim = c(0,max(apply(b1$grad[,], 1, function(t) sum(abs(t))))), ylab = "Complete gradient")
lines(apply(b1$grad[index,], 1, function(t) sum(abs(t))), col = "red")
lines(apply(b10$grad[index,], 1, function(t) sum(abs(t))), col = "blue")
lines(apply(b100$grad[index,], 1, function(t) sum(abs(t))), col = "orange")
legend("topright", col = c("black", "red", "blue", "orange"), legend = c("GA", "SGA", "MB-SGA-10", "MB-SGA-100"), lty = 1)


ga <- SGA(y, X, batch.size = n, t = .05, max.iter = 1e3, adapt = TRUE)
b1 <- SGA(y, X, batch.size = 1, t = 1, max.iter = ga$iter, adapt = TRUE) 
b10 <- SGA(y, X, batch.size = 10, t = 1, max.iter = ga$iter, adapt = TRUE)
b100 <- SGA(y, X, batch.size = 100, t = 1, max.iter = ga$iter, adapt = TRUE) 

index <- 1:1000
plot(apply(ga$grad[index,], 1, function(t) sum(abs(t))), type = 'l', ylim = c(0,max(apply(b1$grad[,], 1, function(t) sum(abs(t))))), ylab = "Complete gradient")
lines(apply(b1$grad[index,], 1, function(t) sum(abs(t))), col = "red")
lines(apply(b10$grad[index,], 1, function(t) sum(abs(t))), col = "blue")
lines(apply(b100$grad[index,], 1, function(t) sum(abs(t))), col = "orange")
legend("topright", col = c("black", "red", "blue", "orange"), legend = c("GA", "SGA", "MB-SGA-10", "MB-SGA-100"), lty = 1)
