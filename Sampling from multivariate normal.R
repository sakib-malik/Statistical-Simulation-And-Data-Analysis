####################################################################
## Generate from multivariate normal distribution
####################################################################
set.seed(1)

par(mfrow = c(2,2))


# Function produces samples from a multivariate normal
multinorm <- function(mu, Sigma, N = 5e2)
{
  decomp <- eigen(Sigma)
  
  # Finding matrix square-root
  Sig.sq <- decomp$vectors %*% diag(decomp$values^(1/2)) %*% solve(decomp$vectors)
  
  samp <- matrix(0, nrow = N, ncol = 2)
  for(i in 1:N)
  {
    Z <- rnorm(2)
    samp[i, ] <- mu + Sig.sq %*%Z
  }
  return(samp)
}


N <- 5e2  # 500 samples

###
# First: Mean (-5, 10) and .5 correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, .5, .5, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)

par(mfrow = c(1,3))
plot(samp, asp = 1, main = "Correlation = .5", xlab = "x_1", ylab = "x_2")
plot(density(samp[,1]), main = "Marginal density for X1")
plot(density(samp[,2]), main = "Marginal density for X2")


par(mfrow = c(2,2))
plot(samp, asp = 1, main = "Correlation = .5", xlab = "x_1", ylab = "x_2")

### 
# Second: Mean (-5, 10) and .99 correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, .99, .99, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)
plot(samp, asp = 1, main = "Correlation = .99", xlab = "x_1", ylab = "x_2")


### 
# Third: Mean (-5, 10) and -.8 correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, -.8, -.8, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)
plot(samp, asp = 1, main = "Correlation = -.8", xlab = "x_1", ylab = "x_2")


### 
# Fourth: Mean (-5, 10) and no correlation
mu <- c(-5,10)
Sigma <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
samp <- multinorm(mu = mu, Sigma = Sigma)
plot(samp, asp = 1, main = "Correlation = 0", xlab = "x_1", ylab = "x_2")


