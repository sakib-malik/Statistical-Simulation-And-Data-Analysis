
#################################
### Plotting weights and
### Weighted average
##################################

set.seed(1)

alpha <- 2
beta <- 5

alpha.p <- alpha + 1
lambda <- 4 #proposal

N <- 1e4
samp <- rgamma(N, shape = alpha.p, rate = lambda)  # importance samples
weights <- dgamma(samp, shape = alpha, rate = beta) /dgamma(samp, shape = alpha.p, rate = lambda)


## Visualizing the target and the importance densities
index <- 1:200
foo <- seq(0, 2, length = 1e3)
plot(foo, dgamma(foo, shape = alpha, rate = beta), 
     type = 'l', col = "black", ylab = "density", ylim = c(0, 3.5), xlab = "Z")
lines(foo, dgamma(foo, shape = alpha.p, rate = lambda), col = "red")
segments(x0 = samp[index], x1 = samp[index], y0 = 0, y1 = weights[index], col = adjustcolor("blue", alpha.f = .2))
legend("topright", legend = c("Target", "Proposal", "Weights"), col = c(1,2,"blue"), lty = 1)





###########################################
### Optimal importance sampling from Gamma
###########################################
set.seed(1)

# Function does importance sampling to estimate second moment of a gamma distribution
imp_gamma <- function(N = 1e3, alpha = 4, beta = 10, moment = 2, imp.alpha = alpha + moment)
{
  fn.value <- numeric(length = N)  
  
  draw <- rgamma(N, shape = imp.alpha, rate = beta) # draw importance samples
  fn.value <- draw^moment * dgamma(draw, shape = alpha, rate = beta) / dgamma(draw, shape = imp.alpha, rate = beta)
  
  return(fn.value)  #return all values
}

N <- 1e4
# Estimate 2nd moment from Gamma(4, 10) using Gamma(4, 10)
# this is IID Monte Carlo
imp_samp <- imp_gamma(N = N, imp.alpha = 4)
mean(imp_samp)
var(imp_samp)

# Estimate 2nd moment from Gamma(4, 10) using Gamma(6, 10)
# this is the optimal proposal
imp_samp <- imp_gamma(N = N)
mean(imp_samp)
var(imp_samp)

# why is the estimate good
foo <- seq(0.001, 5, length = 1e3)
plot(foo, dgamma(foo, shape = 4, rate = 10), type= 'l', ylab = "Density")
lines(foo, dgamma(foo, shape = 6, rate = 10), col = "red")
legend("topright", col = 1:2, lty = 1, legend = c("Reference", "Optimal"))



# Choosing a horrible proposal
# Estimate 2nd moment from Gamma(4, 10) using Gamma(100, 10)
imp_samp <- imp_gamma(N = N, imp.alpha = 100)
mean(imp_samp)   ## estimate is bad too
var(imp_samp)

# why is the estimate bad?
foo <- seq(0.001, 17, length = 1e3)
plot(foo, dgamma(foo, shape = 4, rate = 10), type= 'l', ylab = "Density")
lines(foo, dgamma(foo, shape = 100, rate = 10), col = "red")
legend("topright", col = 1:2, lty = 1, legend = c("Reference", "Importance"))






# This part is not in the notes
# Doing a simulation study for this
# Repeat the above many times to estimate the variability in the estimators
reps <- 1e3
N <- 1e4

var_ests <- matrix(0, nrow = reps, ncol = 4)
mean_ests <- matrix(0, nrow = reps, ncol = 4)
colnames(var_ests) <- c("4", "6", "7", "100")
colnames(mean_ests) <- c("4", "6", "7", "100")

for(i in 1:reps)
{
  imp_samp <- imp_gamma(N = N, imp.alpha = 4)
  mean_ests[i, 1] <- mean(imp_samp)
  var_ests[i, 1] <- var(imp_samp)
  
  imp_samp <- imp_gamma(N = N, imp.alpha = 6)
  mean_ests[i, 2] <- mean(imp_samp)
  var_ests[i, 2] <- var(imp_samp)
  
  imp_samp <- imp_gamma(N = N, imp.alpha = 7)
  mean_ests[i, 3] <- mean(imp_samp)
  var_ests[i, 3] <- var(imp_samp)
  
  imp_samp <- imp_gamma(N = N, imp.alpha = 100)
  mean_ests[i, 4] <- mean(imp_samp)
  var_ests[i, 4] <- var(imp_samp)
}

colMeans(mean_ests)  # Last estimate is horrible
colMeans(var_ests)   # Smallest for 6

