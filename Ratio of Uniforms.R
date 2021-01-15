##################################################
### Ratio of Uniforms for Exp(1)
##################################################
set.seed(1)
# function to sample from the rectangle
drawFromRect <- function(a, b, c)
{
  u <- runif(1, min = 0, max = a)
  v <- runif(1, min = b, max = c)
  return(c(u,v))
}
# sqrt f function
sqrt.f <- function(x) exp(-x/2)

# Starting the process for Exp(1)
a <- 1
b <- 0
c <- 2*exp(-1)
prob.of.acceptance <- 1/(2*a*(c-b))  # true prob. of acceptance for AR

N <- 1e4 # number of samples
samp <- numeric(length = N)
i <- 1
counter <- 0  # to check acceptance
while(i <= N)
{
  counter <- counter + 1
  prop <- drawFromRect(a = a, b = b, c = c)
  vbyu <- prop[2]/prop[1]
  if( prop[1] < sqrt.f(vbyu))
  {
    samp[i] <- vbyu
    i <- i + 1
  }
}

plot(density(samp), main = "Estimated density for Exp(1)")
lines(density(rexp(1e4, 1)), col = "red")
legend("topright", col = c("black", "red"), lty = 1, legend = c("RoU", "Truth"))

(prob.of.acceptance)
# [1] 0.6795705

N/counter  # very close
# [1] 0.6796248


