################################################
## EM Algorithm for the Censored exponential data
################################################

# First we will simulate data
set.seed(10)
n <- 100
lam <- .45
z <- rexp(n, rate = lam)
T <- 1  

# Observed
obs.z <- z[(z > T)]
m <- n - length(obs.z) # to find m
hist(z, breaks = 30, main = "Complete failure times with observed times as dots")  
points(obs.z,rep(0, length(obs.z)), col = 2, pch = 16)  
abline(v = 1, col = 2)



# EM Algorithm
diff <- 100
tol <- 1e-5
iter <- 0
current <- 1/mean(obs.z)  # good starting values
lam.k <- current
store.z <- T
while(diff > tol)
{
  iter <- iter + 1
  # E-step
   Estep <- 1/current - T*exp(-current * T)/(1 - exp(-current * T))
   
   # Storing Zs
  store.z <- c(store.z, Estep)
  
  #Update
  update <- n/(sum(obs.z) + m*Estep)
  
  diff <- abs(current - update)
  current <- update
  lam.k <- c(lam.k, update)
}

current
(false.mle <- (n-m)/sum(obs.z))  # MLE if you ignore the Es (bad MLE)

# Estimate of Z_i not observed
(Z_unobs <- 1/current - T*exp(-current * T)/(1 - exp(-current * T)))



################################################
## Monte Carlo EM Algorithm for the 
## Censored exponential data
################################################
set.seed(1)
#function draws from truncated exponential
trexp <- function(nsim, T, lambda)
{
  count <- 1
  samp <- numeric(length = nsim)
  while(count < nsim)
  {
    draw <- rexp(1, rate = lambda)
    if(draw <= T)
    {
      samp[count] <-  draw
      count <- count + 1
    } else{
      next
    }
  }
  return(samp)
}

diff <- 100
tol <- 1e-5
iter <- 0
current <- 2
lam.k <- current
store.z <- T
mc.size <- 100  # Monte Carlo size
while(diff > tol)
{
  iter <- iter + 1
  
  # E-step with Monte Carlo
  
  trunExp.samples <- trexp(nsim = mc.size, T = T, lambda = current)
  Estep <- mean(trunExp.samples)
  #Estep <- 1/current - T*exp(-current * T)/(1 - exp(-current * T))
  
  # Storing Zs
  store.z <- c(store.z, Estep)
  
  #Update
  update <- n/(sum(obs.z) + m*Estep)
  
  diff <- abs(current - update)
  current <- update
  lam.k <- c(lam.k, update)
}

current




################################################
## Monte CarloEM Algorithm for the 
## Censored exponential data with larger MC size
################################################

diff <- 100
tol <- 1e-5
iter <- 0
current <- 2
lam.k <- current
store.z <- T
mc.size <- 1000  # Monte Carlo size
while(diff > tol)
{
  iter <- iter + 1
  
  # E-step with Monte Carlo
  
  trunExp.samples <- trexp(nsim = mc.size, T = T, lambda = current)
  Estep <- mean(trunExp.samples)
  #Estep <- 1/current - T*exp(-current * T)/(1 - exp(-current * T))
  
  # Storing Zs
  store.z <- c(store.z, Estep)
  
  #Update
  update <- n/(sum(obs.z) + m*Estep)
  
  diff <- abs(current - update)
  current <- update
  lam.k <- c(lam.k, update)
}

current
