#########################################
## Coin Probability
## Posterior distribution
#########################################
set.seed(1)
# generate data
p <- .3

n <- 10
y <- rbinom(n = n, size = 1, prob = p)
x.ax = seq(0, 1, 0.01)
# plotting prior, and posterior

plot(x.ax, rep(p, length(x.ax)), type = 'l', col = "red", ylim = c(0,8), xlab = "p", ylab = "Density", main = "Uniform prior")
lines(x.ax, dbeta(x.ax, shape1 = sum(y) + 1, shape2 = n - sum(y) + 1), col = "blue", lty = 2)
legend("topright", legend = c("Prior", "Posterior for n = 10", "Posterior for n = 100"), col = c("red", "blue", "blue"), lty = c(1,2,2))

# increasing sample size
n <- 100
y <- rbinom(n = n, size = 1, prob = p)

# posterior at larger sample size
lines(x.ax, dbeta(x.ax, shape1 = sum(y) + 1, shape2 = n - sum(y) + 1), col = "blue", lty = 3)



# Changing to Beta prior centered at 1/2
# plotting prior, and posterior
n <- 10
y <- rbinom(n = n, size = 1, prob = p)

plot(x.ax, dbeta(x.ax, shape1 =  5, shape2 = 5), type = 'l', col = "red", ylim = c(0,8), xlab = p, ylab = "Density", main = "Beta(5,5) prior")
lines(x.ax, dbeta(x.ax, shape1 = sum(y) + 5, shape2 = n - sum(y) + 5), col = "blue", lty = 2)

# increasing sample size
n <- 100
y <- rbinom(n = n, size = 1, prob = p)

# posterior at larger sample size
lines(x.ax, dbeta(x.ax, shape1 = sum(y) + 5, shape2 = n - sum(y) + 5), col = "blue", lty = 3)
legend("topright", legend = c("Prior", "Posterior for n = 10", "Posterior for n = 100"), col = c("red", "blue", "blue"), lty = c(1,2,2))


plot(x.ax, dbeta(x.ax, shape1 = sum(y) + 5, shape2 = n - sum(y) + 5), col = "blue", lty = 2, ylim = c(0,8), xlab = p, ylab = "Density", main = "Beta(5,5) prior", type = 'l')
abline(v = qbeta(c(.025, .975), shape1 = sum(y) + 5, shape2 = n - sum(y) + 5), col = "blue", lty = 1)
abline(v = .5, col = "red")
legend("topright", legend = c("Posterior", "95% Credible Interval", "Fair coin, p = .5"), col = c("blue", "blue", "red"), lty = c(2,1,1))


