###########################################
## Linchpin variable sampler
## for Bayesian linear regression for cars
###########################################
set.seed(1)

# loading the dataset
data(cars)
n <- dim(cars)[1]
X <- cbind(1, cars$speed)
y <- cars$dist
p <- dim(X)[2]
a <- 1 # prior parameters
b <- 1 # prior parameters


# We implement Monte Carlo sampling using linchping
N <- 1e4
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


par(mfrow = c(1,3))
plot(density(sig2), main = expression(sigma^2))
plot(density(beta[,1] ), main = expression(beta[1]))
plot(density(beta[,2] ), main = expression(beta[2]))

poster <- cbind(sig2, beta)
colMeans(poster)
apply(poster, 2, quantile, c(.025, .975))


