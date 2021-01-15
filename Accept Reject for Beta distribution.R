#########################
### Accept-reject for
## Beta(4,3) distribution
#########################
set.seed(1)
beta_ar <- function() 
{
	c <- 60 *(3/5)^3 * (2/5)^2
	accept <- 0
	counter <- 0   # count the number of loop
	while(accept == 0)
	{
		counter <- counter + 1
		U <- runif(1)
		prop <- runif(1)

		ratio <- dbeta(prop, shape1 = 4, shape2 = 3)/c

		if(U <= ratio)
		{
			accept <- 1
			return(c(prop, counter))
		}
	}
}

N <- 1e4
samp <- numeric(length = N)
counts <- numeric(length = N)
for(i in 1:N)
{
	foo <- beta_ar()  # I use foo as a dummy name
	samp[i] <- foo[1]
	counts[i] <- foo[2]
}

x <- seq(0, 1, length = 500)
plot(density(samp), main = "Estimated density from 1e4 samples")
lines(x, dbeta(x, shape1 = 4, shape2 = 3), col = "red", lty = 2)
legend("topleft", lty = 1:2, col = c("black", "red"), legend = c("AR", "truth"))

# This is c
(c <- 60 *(3/5)^3 * (2/5)^2)

# This is the mean number of loops required
mean(counts)

#They are almost the same!




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

4/pi
mean(counts)  # very close


