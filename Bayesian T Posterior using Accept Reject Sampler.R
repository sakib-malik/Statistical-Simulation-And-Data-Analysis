###########################################
## Accept-reject to sample from posterior distrbituon
# from t likelihood and normal priors
# Code takes some time to run
###########################################
library(MASS)
set.seed(10)

# Log posterior
log.post <- function(y, mu, nu)
{
  -mu^2/2 - (nu + 1)/2 *sum( log( 1 + (y-mu)^2/nu  ) )
}

# Accept-reject for a given bound M
AR_tmodel <- function(N = 1e2, M = 1)
{
  count <- 0
  attempts <- 0
  samp <- numeric(length = N)
  while(count < N)
  {
    attempts <- attempts + 1
    prop <- rnorm(1)
    
    ratio <- log.post(y = y, mu = prop, nu = nu) + prop^2/2 - log(M)
    if(exp(ratio) > 1) print(exp(ratio))  # Making sure M is correct
    if(runif(1) < exp(ratio))
    {
      count <- count + 1
      samp[count] <- prop
      if(count%% (N) == 0) print(paste("Accepted = ",count, ", Accept Prob. = ", count/attempts))
    }
  }
  return(samp)
}

# Generate the data set.
# Small n and mu close to zero (the prior)
nu <- 3
n <- 10
mu <- 0
y <- mu + rt(n, df = nu)   # Generate the data

mle <- fitdistr(y, "t", df = 3)$estimate[1]    # Find the MLE to construct the tighter upper bound for pi(x)/g(x)
M2 <- prod( (1 + (y - mle)^2/nu ) )^(-(nu+1)/2) + 1e-5 # adding a little to remove numerical approximation errors


# Run the A-R sampler using the two different upper bounds
system.time(out1 <- AR_tmodel(N = 1e4, M = 1) )  # About 50 seconds
system.time(out2 <- AR_tmodel(N = 1e4, M = M2) )  # About .11 seconds

# Posterior mean estimates
c(mean(out1), mean(out2))

# 95% credible interval
quantile(out1, c(.025, .975))
quantile(out2, c(.025, .975))

# Compare the performance
plot(density(rnorm(1e5)), col = "red", type = 'l', ylim = c(0,1.1), main = "Prior and Posterior")#prior
lines(density(out1)) # samples posterior 2
lines(density(out2), col = "blue") # sampled posterior 2
legend("topright", legend = c("Prior", "Posterior 1", "Posterior 2"), col = c("red", "black", "blue"), lty = 1)




# Increase n
n <- 20
y <- mu + rt(n, df = nu)   # Generate the data

mle <- fitdistr(y, "t", df = 3)$estimate[1]    # Find the MLE to construct the tighter upper bound for pi(x)/g(x)
M2 <- prod( (1 + (y - mle)^2/nu ) )^(-(nu+1)/2) + 1e-6


# Run the A-R sampler using the two different upper bounds
system.time(out1 <- AR_tmodel(N = 1e1, M = 1) )  # too slow for even 10 samples
system.time(out2 <- AR_tmodel(N = 1e4, M = M2) )  # About .11 seconds

# Posterior mean estimates
mean(out2)

# 95% credible interval
quantile(out2, c(.025, .975))

# Compare the performance. Not enough samples
plot(density(rnorm(1e5)), col = "red", type = 'l', ylim = c(0,1))#prior
lines(density(out2), col = "blue") # sampled posterior






# Repeating the above experiment with true mean far away from 0
# Generate the data set.
# Small n and true mu farther from zero, the prior mean
nu <- 3
n <- 5
mu <- 5
y <- mu + rt(n, df = nu)   # Generate the data

mle <- fitdistr(y, "t", df = 3)$estimate[1]    # Find the MLE to construct the tighter upper bound for pi(x)/g(x)
M2 <- prod( (1 + (y - mle)^2/nu ) )^(-(nu+1)/2) + 1e-4


# Run the A-R sampler using the two different upper bounds
system.time(out1 <- AR_tmodel(N = 5, M = 1) )  # About 6 seconds
system.time(out2 <- AR_tmodel(N = 5, M = M) )  # About 5 seconds

