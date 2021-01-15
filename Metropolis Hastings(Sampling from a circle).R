###############################################
## Accept-reject for uniform distribution
## on a circle
###############################################
set.seed(1)
circle_ar <- function() 
{
  accept <- 0
  counter <- 0   # count the number of loop
  while(accept == 0)
  {
    counter <- counter + 1
    prop <- runif(2, min = -1, max = 1) #c(U1, U2)
    
    in.or.out <- prop[1]^2 + prop[2]^2 < 1
    
    if(in.or.out)
    {
      accept <- 1
      return(c(prop, counter))
    }
  }
}

N <- 1e4
samp <- matrix(0, ncol = 2, nrow = N)
counts <- numeric(length = N)
for(i in 1:N)
{
  foo <- circle_ar()  # I use foo as a dummy name
  samp[i,] <- foo[1:2]
  counts[i] <- foo[3]
}

plot(samp[,1], samp[,2], xlab = "x", ylab = "y",main = "Uniform samples from a circle", asp = 1)



###########################################
## MCMC samples from a uniform distribution
## within a circle
###########################################

N <- 1e4
xsamp <- numeric(length = N)
ysamp <- numeric(length = N)
xprop <- numeric(length = N)
yprop <- numeric(length = N)
acc <- numeric(length = N)
acc[1] <- 1
samp_eval <- numeric(length = N)
colors <- ifelse(samp_eval, "blue", "red")
  
h <- .2   #change h and see what happens
xsamp[1] <- .2
xprop[1] <- .2
ysamp[1] <- .6
yprop[1] <- .6
for(i in 2:N)
{
  propx <- xsamp[i-1] + runif(1, -h, h)
  propy <- ysamp[i-1] + runif(1, -h, h)
  
  acc[i] <- (propx^2 + propy^2) < 1
  
  if(acc[i] == 0)
  {
    xsamp[i] <- xsamp[i-1]
    ysamp[i] <- ysamp[i-1]
  } else
  {
    xsamp[i] <- propx
    ysamp[i] <- propy
  }
  
  xprop[i] <- propx
  yprop[i] <- propy
}
ind <- xsamp > 0 

samp_eval <- (xsamp)^2 + (ysamp)^2 < 1
colors <- ifelse(acc, "blue", "red")

x <- seq(-1.1,1.1, length = 1e4)
y <- seq(-1.1,1.1, length = 1e4)

for(i in 1:100)
{
  plot(x, y, type = "n", col = "red", main = "MCMC sampling from a circle", ylim = range(y), xlim = range(x), asp = 1)
  draw.circle(0,0,radius = 1)
  points(xsamp[i], ysamp[i], col = "green", pch = 16, cex = 1) 
  segments(x0 = xsamp[i]-h, y0 = ysamp[i]-h, x1 = xsamp[i]+h, y1 = ysamp[i]-h, lwd = 2)
  segments(x0 = xsamp[i]-h, y0 = ysamp[i]-h, x1 = xsamp[i]-h, y1 = ysamp[i]+h, lwd = 2)
  segments(x0 = xsamp[i]-h, y0 = ysamp[i]+h, x1 = xsamp[i]+h, y1 = ysamp[i]+h, lwd = 2)
  segments(x0 = xsamp[i]+h, y0 = ysamp[i]+h, x1 = xsamp[i]+h, y1 = ysamp[i]-h, lwd = 2)
  points(xsamp[1:i-1], ysamp[1:i-1], col = alpha("blue", .3), pch = 16, cex = 1) 
  points(xprop[i+1], yprop[i+1], col = colors[i+1], pch = 16, cex = 1) 
  
  Sys.sleep(.20)
}

plot(x, y, type = "n", col = "red", main = "MCMC sampling from a circle", ylim = range(y), xlim = range(x), asp = 1)
draw.circle(0,0,radius = 1)
points(xprop, yprop, col = colors, pch = 16, cex = 1) 

# Trace plot and ACF plot
par(mfrow = c(1,2))
acf(xsamp, main = paste("Autocorrelation for x for h = ", h))
plot(xsamp, type = 'l', main = paste("Trace plot for x for h = ", h))


