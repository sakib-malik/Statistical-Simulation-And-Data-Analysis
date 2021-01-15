###########################################
## MCMC for Bayesian linear regression
## for the cars dataset
###########################################
set.seed(1)

data(cars)
n <- dim(cars)[1]
X <- cbind(1, cars$speed)
y <- cars$dist

# Calculates the log-posterior since if not on the log scale, numerical instability
log.post <- function(y, X, beta, s2, a, b, mu = 0)
{
  n <- length(y)
  p <- dim(X)[2]
  foo <- (-n/2 - p/2 - a -1) * log(s2) + (-sum((y - X%*%beta)^2)/2 - sum((beta - mu)^2)/2 - b)/s2 
  return(foo)
}

# MCMC function with three different step-sizes for each component
BLRmcmc <- function(y, X,T = 1e4, a, b, mu = 0, start, h1 = 3, h2 = 3, h3 = 3)
{
  p <- dim(X)[2]

  beta <- matrix(0, nrow = T, ncol = p)
  s2 <- numeric(length = T)
  
  beta[1,] <- start[1:p]   # Set starting values
  s2[1] <- start[p+1]
  
  acc <- 0
  for(t in 2:T)
  {
    prop <- c(beta[t-1,], s2[t-1]) + rnorm(p+1, sd = sqrt(c(h1,h2,h3)))   # Proposal
    print(prop)
    if(prop[p+1] < 0)   # If the variance is negative, it's a straight reject
    {
      beta[t,] <- beta[t-1,]
      s2[t] <- s2[t-1]
    } else{
      
      # difference of log posteriors
      ratio <- log.post(y = y, X = X, beta = prop[1:p], s2 = prop[p+1], a = a, b = b, mu = mu) - log.post(y = y, X = X, beta = beta[t-1,], s2 = s2[t-1], a = a, b = b, mu = mu)
      if(runif(1) < exp(ratio))
      {
        acc <- acc + 1
        beta[t,] <- prop[1:p]
        s2[t] <- prop[p+1]
      } else{
       beta[t, ] <- beta[t-1, ]
       s2[t] <- s2[t-1]
      }
    }
  }
  print(paste("Acceptance prob = ", acc/T))
  return(list("beta" = beta, "s2" = s2))
}


lm.obj <- lm(y ~ X - 1)
s2.start <- var(lm.obj $residuals)

# Three different step sizes for three different components. Change h1, h2, h3 to see how Acceptance prob is affected
# Tuned to be aroun 32% since dimension is 3
# Also notice the starting values
chain <- BLRmcmc(T = 1e5, y = y, X = X, a = 1, b = 1, start = c(lm.obj$coeff, s2.start), h1 = 20, h2 = .35, h3 = 100)

# Poster mean Estimates and credible intervals
post.mean <- c(colMeans(chain$beta), mean(chain$s2))
post.quant1 <- c(quantile(chain$beta[,1], .05), quantile(chain$beta[,2], .05), quantile(chain$s2, .05))
post.quant2 <- c(quantile(chain$beta[,1], .95), quantile(chain$beta[,2], .95), quantile(chain$s2, .95))


# Estimates of intercept
c(post.quant1[1],post.mean[1], post.quant2[1])
lm.obj$coefficients[1]  # Compare with MLE

# Estimates of coefficient of speed
c(post.quant1[2], post.mean[2], post.quant2[2])
lm.obj$coefficients[2]  # Compare with MLE

# Estimates of variance
c(post.quant1[3],post.mean[3],  post.quant2[3])
var(lm.obj $residuals)  # Compare with MLE


#Visualizing plots
par(mfrow = c(3,3))
plot.ts(chain$beta[,1])
plot.ts(chain$beta[,2])
plot.ts(chain$s2)

acf(chain$beta[ ,1])
acf(chain$beta[ ,2])
acf(chain$s2)

plot(density(chain$beta[,1]))
abline(v = c(post.mean[1], post.quant1[1], post.quant2[1]), col = c("red", "blue", "blue"))
plot(density(chain$beta[,2]))
abline(v = c(post.mean[2], post.quant1[2], post.quant2[2]), col = c("red", "blue", "blue"))
plot(density(chain$s2))
abline(v = c(post.mean[3], post.quant1[3], post.quant2[3]), col = c("red", "blue", "blue"))






