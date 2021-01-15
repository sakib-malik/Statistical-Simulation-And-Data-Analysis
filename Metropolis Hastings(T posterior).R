###########################################
## MH samples from posterior distrbituon
# from t likelihood and normal prior
###########################################
set.seed(10)
log.post <- function(y, mu, nu)
{
  -mu^2/2 - (nu + 1)/2 *sum( log( 1 + (y-mu)^2/nu  ) )
}

# Generate the data set.
# No need for small n and mu close to zero (the prior)
nu <- 3
n <- 50
mu <- 4
y <- mu + rt(n, df = nu)   # Generate the data


# Now compare to MCMC for the same target distribution
T <- 1e5
mc.samp <- numeric(length = T)
mc.samp[1] <- mean(y)
acc <- 0
foo <- proc.time()
for(t in 2:T)
{
  prop <- rnorm(1, mean = mc.samp[t-1], sd = .4)
  ratio <- log.post(y = y, mu = prop, nu = nu) - log.post(y = y, mu = mc.samp[t-1], nu = nu)
  
  if(runif(1) < exp(ratio))
  {
    mc.samp[t] <- prop
    acc <- acc + 1
  }else{
    mc.samp[t] <- mc.samp[t-1]
  }
}
proc.time() - foo
print(acc/T)

par(mfrow = c(1,3))
plot(density(mc.samp), col = "blue", xlim = c(-3,5), main = "Density plot")
lines(density(rnorm(1e5)), col = "red")
legend("topleft", col = c("red", "blue"), legend = c("Prior", "Posterior"), lty = 1)
acf(mc.samp, main = "ACF Plot")
plot.ts(mc.samp, main = "Trace plot")

# Since the sampler looks good, we can now
# estimate the posterior mean and quantiles
mean(mc.samp)
quantile(mc.samp, c(.025, .975))

# a short run again to indicate what happens
# with a bad starting value
# But choose a BAD starting values
T <- 1e3
mc.samp <- numeric(length = T)
mc.samp[1] <- 50 ## bad starting value
acc <- 0
foo <- proc.time()
for(t in 2:T)
{
  prop <- rnorm(1, mean = mc.samp[t-1], sd = .5)
  ratio <- log.post(y = y, mu = prop, nu = nu) - log.post(y = y, mu = mc.samp[t-1], nu = nu)
  
  if(runif(1) < exp(ratio))
  {
    mc.samp[t] <- prop
    acc <- acc + 1
  }else{
    mc.samp[t] <- mc.samp[t-1]
  }
}

par(mfrow = c(1,1))
plot.ts(mc.samp, main = "Trace plot with bad starting value")


