#####################################
## Simulated Annealing
## Demonstrative example
#####################################
set.seed(1)
fn <- function(x, T = 1)
{
  h <- ( cos(50*x) + sin(20*x) )^2 
  exp(h/T) * (0 < x & x < 1)
}

x <- seq(0, 1, length = 5e2)
plot(x, fn(x), type = 'l', ylim = c(0,150), ylab = "exp(f/T)")
for(t in 2:4)
{
  lines(x, fn(x, T = 1/(log(t))), col = (t))
}
legend("topright", col = 1:4, lty = 1, legend = c("T = 1", "T = .83", "T = .75", "T = .71"))

simAn <- function(N = 10, r = .3)
{
  x <- numeric(length = N)
  x[1] <- runif(1)
  
  for(k in 2:N)
  {
    # U(x - r, x + r)
    a <- runif(1, x[k-1] - r,  x[k-1] + r)
    T <- 1/(log(k))
    
    ratio <- fn(a,T)/fn(x[k-1], T)
    if( runif(1) < ratio)
    {
      x[k] <- a # accept
    } else{
      x[k] <- x[k-1] # reject, so stay
    }
  }
  x
 return(x) 
}

N <- 500
sim <- simAn(N = N)
sim[which.max(fn(sim))]  # theta^*

plot(x, fn(x), type = 'l', ylab = "exp(f/T)")
points(sim, fn(sim), pch = 16, col = 1)


#####################################
## Location Cauchy likelihood example
#####################################

################################################
## MLE for location Cauchy distribution using
## Newton-Raphson method
## We will plot the likelihood as well
################################################
set.seed(1)
mu.star <- 5  # Setting true mu
n <- 4  # sample size
X <- rt(n, df = 1) + mu.star

## Function calculates the exp(like/T)
log.like <- function(mu, X, T = 1)
{
  n <- length(X)
  rtn <- -n*log(pi) - sum( log(1 + (X - mu)^2) )
  return(exp(rtn/T))
}

mu.x <- seq(-10, 40, length = 1e3)  # A sequence of mu's 
ll.est <- log(sapply(mu.x, log.like, X))  # evaluating log-likelihood at the mus
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))  # plotting log-likelihood. Not concave, so we need to choose good starting values.

# Simulated annealing algorithm
simAn <- function(N = 10, r = .5)
{
  x <- numeric(length = N)
  x[1] <- runif(1, min = -10, max = 40)
  fn.value <- numeric(length = N)
  
  fn.value[1] <- log.like(mu = x[1], X, T = 1)
  for(k in 2:N)
  {
    a <- runif(1, x[k-1] - r, x[k-1] + r)
    T <- 1/(1 + log(log(k)))
    ratio <- log.like(mu = a, X, T)/log.like(mu = x[k-1], X, T)
    if( runif(1) < ratio)
    {
      x[k] <- a
    } else{
      x[k] <- x[k-1]
    }
    fn.value[k] <- log.like(mu = x[k], X, T = 1)
  }
  x
  return(list("x" = x, "fn.value" = fn.value)) 
}

par(mfrow = c(2,2))

# Four different runs all converge.
sim <- simAn(N = 1e2, r = 5)
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu)) 
points(sim$x, log(sim$fn.value), pch = 16, col = adjustcolor("blue", alpha = .4))

sim <- simAn(N = 1e2, r = 5)
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))
points(sim$x, log(sim$fn.value), pch = 16, col = adjustcolor("darkred", alpha = .4))

sim <- simAn(N = 1e2, r = 5)
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))
points(sim$x, log(sim$fn.value), pch = 16, col = adjustcolor("darkgreen", alpha = .4))

sim <- simAn(N = 1e2, r = 5)
plot(mu.x, ll.est, type = 'l', ylab = "log-likelihood", xlab = expression(mu))
points(sim$x, log(sim$fn.value), pch = 16, col = adjustcolor("purple", alpha = .4))


par(mfrow = c(1,2))
## Different values of r
# very large r
sim <- simAn(N = 1e3, r = 500)
plot(mu.x, ll.est, type = 'l', main = "r = 500. Many rejections",  ylab = "log-likelihood", xlab = expression(mu))  
points(sim$x, log(sim$fn.value), pch = 16, col = adjustcolor("blue", alpha = .2))

#very small r
plot(mu.x, ll.est, type = 'l', main = "r = .1. Many small acceptances", ylab = "log-likelihood", xlab = expression(mu)) 
sim <- simAn(N = 1e3, r = .1)
points(sim$x, log(sim$fn.value), pch = 16, col = adjustcolor("blue", alpha = .2))
