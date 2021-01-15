###########################################
## Bayesian logistic regression
## with MH implementation
###########################################
# log posterior
logf <- function(beta)
{
  one.minus.yx <- (1 - y)*X	
  -sum(beta^2)/2 - sum(log(1 + exp(-X%*%beta))) - sum(one.minus.yx%*%beta)
}

bayes_logit_mh <- function(y, X, N = 1e4, prop.sd = .35)
{
  p <- dim(X)[2]
  one.minus.yx <- (1 - y)*X	
  
  # starting value is the MLE
  ##foo <- glm(y ~X - 1, family = binomial("logit"))$coef
  foo <- rep(0, p)
  foo <- as.numeric(foo)
  beta <- as.matrix(foo, ncol = 1)
  beta.mat <- matrix(0, nrow = N, ncol = p)
  beta.mat[1, ] <- as.numeric(beta)
  accept <- 0
  
  for(i in 2:N)
  {
    #symmetric density
    prop <- rnorm(p, mean = beta, sd = prop.sd)
    
    # log of the MH ratio
    log.rat <- logf(prop) - logf(beta)
    if(log(runif(1)) < log.rat)
    {
      beta <- prop
      accept <- accept + 1
    }
    beta.mat[i, ] <- beta
  }
  print(paste("Acceptance Prob = ", accept/N))
  return(beta.mat)
}

titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")

y <- titanic[,1]
X <- as.matrix(titanic[, -1])

## acceptance is too low! we want 23%
## so decrease proposal variance
chain <- bayes_logit_mh(y = y, X = X, N = 1e3, prop.sd = .35) 

# still too low
chain <- bayes_logit_mh(y = y, X = X, N = 1e3, prop.sd = .1) 

# now its better
chain <- bayes_logit_mh(y = y, X = X, N = 1e3, prop.sd = .0065) 

# will now run the chain much longer for 10^5
# takes a few seconds
chain <- bayes_logit_mh(y = y, X = X, N = 1e5, prop.sd = .0065) 

# all trace plots
plot.ts(chain, main = "Trace plots")

par(mfrow = c(2,3))
# all ACF plots
for(i in 1:dim(chain)[2])
{
  acf(chain[,i], main = paste("ACF of Comp ", i))
}

# all density plots plots
for(i in 1:dim(chain)[2])
{
  plot(density(chain[,i]), main = paste("Density of Comp ", i))
}

# we see above that some components are ok, but 4 components are
# moving very slowly. This is because we are using the same proposal
# variance for each component, which is not adequate here.
# Below now I use different proposal variances for different
#components.


chain <- bayes_logit_mh(y = y, X = X, N = 1e5, prop.sd = c(.08, .08, .0065, .03, .03, .0065)) 

# all trace plots
plot.ts(chain, main = "Trace plots")

par(mfrow = c(2,3))
# all ACF plots
for(i in 1:dim(chain)[2])
{
  acf(chain[,i], main = paste("ACF of Comp ", i))
}

# all density plots plots
for(i in 1:dim(chain)[2])
{
  plot(density(chain[,i]), main = paste("Density of Comp ", i))
}
for(i in 1:dim(chain)[2]){
  print(paste("Beta",i," : "))
  print(mean(chain[,i]))
  print(quantile(chain[,i], c(0.025, 0.975)))
}
chain[1 ,]
