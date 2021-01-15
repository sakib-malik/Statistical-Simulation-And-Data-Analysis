###########################################
### Importance sampline from Gamma
###########################################
set.seed(1)

alpha <- 2
beta <- 5

k <- 2  # second moment
(truth <- (alpha / beta^2) + (alpha/beta)^2)  # true second moment

lambda <- 3 #proposal

N <- 1e4
samp <- rexp(N, rate = lambda)  # importance samples
func <- samp^k * dgamma(samp, shape = alpha, rate = beta) / dexp(samp, rate = lambda)
mean(func)  # truth is .24


## Visualizing the target and the importance densities
foo <- seq(0, 2, length = 1e3)
plot(foo, dgamma(foo, shape = alpha, rate = beta), 
  type = 'l', col = "black", ylab = "density", ylim = c(0, 3.5))
lines(foo, dexp(foo, rate = lambda), col = "red")
points(x = samp[1:100],y = rep(0, 100), col = "blue")
legend("topright", legend = c("Gamma(2,5)", "Exp(3)"), col = c(1,2), lty = 1)


## Checking convergence
N <- 1e5 # very large N
samp <- rexp(N, rate = lambda)  # importance samples
func <- samp^k * dgamma(samp, shape = alpha, rate = beta) / dexp(samp, rate = lambda)

# Plotting the running average
plot(1:N, cumsum(func)/(1:N), type = 'l', xlab = "N", ylab = "Running average")
abline(h = truth, col = "red")



## Checking if unbiased or not
N <- 1e4
r <- 1e3
ests <- numeric(length = r)
for(a in 1:r)
{
  samp <- rexp(N, rate = lambda)  # importance samples
  func <- samp^k * dgamma(samp, shape = alpha, rate = beta) / dexp(samp, rate = lambda)
  ests[a] <- mean(func)
}
mean(ests - truth)  # very close to 0

# looking at variance
var(ests) # This is var(theta_g) = sigma^2_g/N
N*var(ests)  # pretty small

## When Accept-reject fails
# If lambda > beta, we know accept-reject fails, let's see what the variance is then
lambda <- 10
N <- 1e4
r <- 1e3
ests <- numeric(length = r)
for(a in 1:r)
{
  samp <- rexp(N, rate = lambda)  # importance samples
  func <- samp^k * dgamma(samp, shape = alpha, rate = beta) / dexp(samp, rate = lambda)
  ests[a] <- mean(func)
}
mean(ests - truth)  # very close to 0

var(ests) # This is var(theta_g) = sigma^2_g/N
N*var(ests)  # Variance is much larger now!



## Checking convergence again
# Convergence is not affected, but it takes MUCH longer
# to get good convergence 
par(mfrow = c(1,2))
N <- 1e5 # very large N
samp <- rexp(N, rate = lambda)  # importance samples
func <- samp^k * dgamma(samp, shape = alpha, rate = beta) / dexp(samp, rate = lambda)

# Plotting the running average
plot(1:N, cumsum(func)/(1:N), type = 'l', xlab = "N", ylab = "Running average")
abline(h = truth, col = "red")


N <- 1e6 # very large N
samp <- rexp(N, rate = lambda)  # importance samples
func <- samp^k * dgamma(samp, shape = alpha, rate = beta) / dexp(samp, rate = lambda)

# Plotting the running average
plot(1:N, cumsum(func)/(1:N), type = 'l', xlab = "N", ylab = "Running average")
abline(h = truth, col = "red")

