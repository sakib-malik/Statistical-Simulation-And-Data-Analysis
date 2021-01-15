################################################
## EM Algorithm for the Old Faithful Geyser data
################################################

data(faithful)
head(faithful)

x <- faithful$eruptions
hist(x, freq = FALSE, breaks = 30, main = "Eruptions")

# (pi_1, mu_1, mu_2, sigma^2_1, sigma^2_2)
theta <- c(.6, 1,5, 1, 1)
current <- theta
diff <- 100
tol <- 1e-5
iter <- 0
store <- current

while(diff > tol)
{
  iter <- iter + 1
  
  # E step: find gamma_{i,c,k} for just c = 1, since for c = 2 is just 1-Ep
  Ep <- current[1]*dnorm(x, current[2], sqrt(current[4]))/
    (current[1]*dnorm(x, current[2], sqrt(current[4])) + (1 - current[1])*dnorm(x, current[3], sqrt(current[5])))
  
  # M-step
  theta[1] <- mean(Ep)
  theta[2] <- sum(Ep*x) / sum(Ep)
  theta[3] <- sum((1-Ep)*x) / sum(1-Ep)
  theta[4] <- sum(Ep*(x - theta[2])^2) / sum(Ep)
  theta[5] <- sum((1-Ep)*(x - theta[3])^2) / sum(1-Ep)
  
  diff <- max( abs(theta - current))
  current <- theta
  store <- rbind(store, theta)
}
current

# Final estimates of the probability
# that each observation is in Class C.
Prob.Z <- current[1]*dnorm(x, current[2], sqrt(current[4]))/
  (current[1]*dnorm(x, current[2], sqrt(current[4])) + (1 - current[1])*dnorm(x, current[3], sqrt(current[5])))

head(round(Prob.Z, 4))

# Make plot of iterative model fits
for(i in 1:dim(store)[1])
{
  test.x <- seq(min(x), max(x), length = 1000)
  test.y <- store[i,1]* dnorm(test.x, mean = store[i,2], sd = sqrt(store[i,4])) + (1-store[i,1]) *dnorm(test.x, mean = store[i,3], sd = sqrt(store[i,5]))
  lines(test.x, test.y, col = rgb(1,0,0, alpha = .5))
}
lines(test.x, test.y, col = rgb(0,0,1, alpha = 1))

# add color
color <- 1*(Ep < .5) + 3*(Ep >= .5)
points(x, rep(0, length(x)), pch = 16, col = color)




################################################
## Doing this again on the Waiting times

x <- faithful$waiting

# (pi_1, mu_1, mu_2, sigma^2_1, sigma^2_2)
theta <- c(.6, 30,100, 1, 1)
current <- theta
diff <- 100
tol <- 1e-5
iter <- 0
store <- current

while(diff > tol)
{
  iter <- iter + 1
  
  # E step: find gamma_{i,c,k} for just c = 1, since for c = 2 is just 1-Ep
  Ep <- current[1]*dnorm(x, current[2], sqrt(current[4]))/
    (current[1]*dnorm(x, current[2], sqrt(current[4])) + (1 - current[1])*dnorm(x, current[3], sqrt(current[5])))
  
  # M-step
  theta[1] <- mean(Ep)
  theta[2] <- sum(Ep*x) / sum(Ep)
  theta[3] <- sum((1-Ep)*x) / sum(1-Ep)
  theta[4] <- sum(Ep*(x - theta[2])^2) / sum(Ep)
  theta[5] <- sum((1-Ep)*(x - theta[3])^2) / sum(1-Ep)
  
  diff <- max( abs(theta - current))
  current <- theta
  store <- rbind(store, theta)
}
current

#hist(x, freq = FALSE, breaks = 30, main = "Wiating times") # a little bit of asymmetry in the first eruptions, so maybe Gaussian isn't the right model
for(i in 1:dim(store)[1])
{
  test.x <- seq(min(x), max(x), length = 1000)
  test.y <- store[i,1]* dnorm(test.x, mean = store[i,2], sd = sqrt(store[i,4])) + (1-store[i,1]) *dnorm(test.x, mean = store[i,3], sd = sqrt(store[i,5]))
  lines(test.x, test.y, col = rgb(1,0,0, alpha = .5))
}
lines(test.x, test.y, col = rgb(0,0,1, alpha = 1))


################################################
# Multivariate mixture model on this
################################################
library(mvtnorm)

x <- as.matrix(faithful)

# (pi_1, mu_1, mu_2, sigma^2_1, sigma^2_2)
mu1 <- c(2.8,75)
mu2 <- c(3.6,58)
sigma1 <- matrix(c(.8, 7, 7, 70), ncol = 2, nrow = 2)
sigma2 <- matrix(c(.8, 7, 7, 70), ncol = 2, nrow = 2)
theta <- c(.6, mu1, mu2, sigma1, sigma2)

current <- theta
diff <- 100
tol <- 1e-10
iter <- 0
store <- current


  
  while(diff > tol)
  {
    iter <- iter + 1
    
    sigma1 <- matrix(current[6:9], ncol = 2, nrow = 2)
    sigma2 <- matrix(current[10:13], ncol = 2, nrow = 2)
    
    # E step: find gamma_{i,c,k} for just c = 1, since for c = 2 is just 1-Ep
    Ep <- current[1] * dmvnorm(x, mean = current[2:3], sigma = sigma1)/
      (current[1]*dmvnorm(x, mean = current[2:3], sigma = sigma1) + (1 - current[1])*dmvnorm(x, mean = current[4:5], sigma = sigma2))
    
    # M-step
    theta[1] <- mean(Ep)
    theta[2:3] <- colSums(Ep*x) / sum(Ep)
    theta[4:5] <- colSums((1-Ep)*x) / sum(1-Ep)
    theta[6:9] <- cov.wt(x, Ep)$cov 
    theta[10:13] <- cov.wt(x, 1-Ep)$cov
    
    diff <- max( abs(theta - current))
    current <- theta
    store <- rbind(store, theta)
  }
  
  allot <- 1*(Ep < .5) + 2*(Ep >= .5)  ## Final allotment of classification
  plot(x, col = allot) # plot allotment
  


