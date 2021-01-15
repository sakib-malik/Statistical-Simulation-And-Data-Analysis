###########################################
## Accept Reject algorithm to draw from
## Binomial(n,p)
###########################################
set.seed(1)

# Function draws one value from Binom(n,p)
# n = number of trials
# p = probability of success
draw_binom <- function(n, p)
{
  accept <- 0
  
  # upper bound calculated in the notes
  x <- 0:n
  all_c <- choose(n,x) * (1-p)^(n - 2*x) * p^(x-1)
  c <- max(all_c) + .001 # final c with slight increase for numerical stability.
  
  
  while(accept == 0)
  {
    U <- runif(1)
    prop <- rgeom(1, prob = p) #draw proposal
    
    ratio <- dbinom(x = prop, size = n, prob = p)/
      (c* dgeom(x = prop, prob = p))
    if(U < ratio)
    {
      accept <- 1
      rtn <- prop
    }
  }
  return(rtn)
}

draw_binom(n = 10, p = .25)


###
# If we want X1, ..., Xn ~ Binom(n.p)
# we need to call the function multiple times

# sample size
N <- 1e3
samp <- numeric(N)
for(t in 1:N)
{
  samp[t] <- draw_binom(n = 10, p = .25)
}
mean(samp) #should be n*p = 2.5
