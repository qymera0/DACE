# Loading packages

library(DiceKriging)

library(DiceEval)

library(lhs)

library(sensitivity)

# 4.1 Introdutory Example 1D

inputs <- c(-1, -0.5, 0, 0.5, 1)

outputs <- c(-9, -5, -1, 9, 11)

theta <- 0.4

sigma <- 5

trend <- c(0, 11, 2)

model <- km(formula = ~x + I(x^2),
            design = data.frame(x = inputs),
            response = outputs,
            covtype = "matern5_2",
            coef.trend = trend,
            coef.cov = theta,
            coef.var = sigma^2)

model

# predict new data

t <- seq(from = -2, to = 2, length = 200)

p <- predict(model, newdata = data.frame(x = t), type = "SK")

# Plot results

plot(t, p$mean, 
     type = "l",
     xlim = c(-2, 2),
     ylim = c(-30, 30),
     xlab = "x",
     ylab = "y")

lines(t, p$lower95, col = "black", lty = 2)

lines(t, p$upper95, col = "black", lty = 2)

points(inputs, outputs, col = "red", pch = 19)

abline(h = 0)

# Influence of Range Parameters

trend <- 0

model <- km(formula = ~1,
            design = data.frame(x = inputs),
            response = outputs,
            covtype = "matern5_2",
            coef.trend = trend,
            coef.cov = theta,
            coef.var = sigma^2)

model

# predict new data

t <- seq(from = -2, to = 2, length = 200)

p <- predict(model, newdata = data.frame(x = t), type = "SK")

# Plot results

plot(t, p$mean, 
     type = "l",
     xlim = c(-2, 2),
     ylim = c(-30, 30),
     xlab = "x",
     ylab = "y")

lines(t, p$lower95, col = "black", lty = 2)

lines(t, p$upper95, col = "black", lty = 2)

points(inputs, outputs, col = "red", pch = 19)

abline(h = 0)

# Afine trend

trend <- c(0, 10)

model <- km(formula = ~x,
            design = data.frame(x = inputs),
            response = outputs,
            covtype = "matern5_2",
            coef.trend = trend,
            coef.cov = theta,
            coef.var = sigma^2)

model

# predict new data

t <- seq(from = -2, to = 2, length = 200)

p <- predict(model, newdata = data.frame(x = t), type = "SK")

# Plot results

plot(t, p$mean, 
     type = "l",
     xlim = c(-2, 2),
     ylim = c(-30, 30),
     xlab = "x",
     ylab = "y")

lines(t, p$lower95, col = "black", lty = 2)

lines(t, p$upper95, col = "black", lty = 2)

points(inputs, outputs, col = "red", pch = 19)

abline(h = 0)

# Affine trend with diferent formula

trend <- c(0, 10)

model <- km(formula = ~.,
            design = data.frame(x = inputs),
            response = outputs,
            covtype = "matern5_2",
            coef.trend = trend,
            coef.cov = theta,
            coef.var = sigma^2)

model

# predict new data

t <- seq(from = -2, to = 2, length = 200)

p <- predict(model, newdata = data.frame(x = t), type = "SK")

# Plot results

plot(t, p$mean, 
     type = "l",
     xlim = c(-2, 2),
     ylim = c(-30, 30),
     xlab = "x",
     ylab = "y")

lines(t, p$lower95, col = "black", lty = 2)

lines(t, p$upper95, col = "black", lty = 2)

points(inputs, outputs, col = "red", pch = 19)

abline(h = 0)

# Sine trend

trend <- c(1, 15)

model <- km(formula = ~sin(pi/4 * x),
            design = data.frame(x = inputs),
            response = outputs,
            covtype = "matern5_2",
            coef.trend = trend,
            coef.cov = theta,
            coef.var = sigma^2)

model

# predict new data

t <- seq(from = -2, to = 2, length = 200)

p <- predict(model, newdata = data.frame(x = t), type = "SK")

# Plot results

plot(t, p$mean, 
     type = "l",
     xlim = c(-2, 2),
     ylim = c(-30, 30),
     xlab = "x",
     ylab = "y")

lines(t, p$lower95, col = "black", lty = 2)

lines(t, p$upper95, col = "black", lty = 2)

points(inputs, outputs, col = "red", pch = 19)

abline(h = 0)

# Simulations

t <- seq(from = -2, to = 2, length = 200)

y <- simulate(model, nsim = 5, newdata = data.frame (x = t)) # each row has one simulation

# plot simulations

trend <- c(0, 11, 2)

ytrend <- trend[1] + trend[2] * t + trend[3] * t^2

par(mfrow = c(1, 1))

plot(t, ytrend,
     type = "l",
     col = "black",
     ylab = "y",
     lty = "dashed",
     ylim = c(min(ytrend) - 2 * sigma, max(ytrend) + 2 * sigma))
 
for(i in 1:5) lines (t, y[i, ], col = i)

# Conditional Simulation

t <- seq(from = -2, to = 2, length = 200)

y <- simulate(model, nsim = 5, newdata = data.frame (x = t), cond = TRUE) # each row has one simulation

# plot simulations

trend <- c(0, 11, 2)

ytrend <- trend[1] + trend[2] * t + trend[3] * t^2

par(mfrow = c(1, 1))

plot(t, ytrend,
     type = "l",
     col = "black",
     ylab = "y",
     lty = "dashed",
     ylim = c(min(ytrend) - 2 * sigma, max(ytrend) + 2 * sigma))

for(i in 1:5) lines (t, y[i, ], col = i)

points(inputs, outputs, col = "red", pch = 19)

# 4.3 Estimation and Validation of Kriging Models

rm(list = ls())

X <- expand.grid(x1 = seq(0, 1, length = 4), x2 = seq(0, 1, length = 4))

y <- apply(X, 1, branin)

m <- km(~ ., 
        design = X, 
        response = y, 
        covtype = "gauss",
        control = list(pop.size = 50,
                       trace = FALSE)) # Increased initial search points to 50

m <- km(~ ., 
        design = X, 
        response = y, 
        covtype = "gauss",
        optim.method = "gen",
        control = list(trace = FALSE))

m
        
# Plot of loglikellyhood function

n.grid <- 30

x.grid <- seq(0.01, 2, length = n.grid)

X.grid <- expand.grid(x.grid, x.grid)

logLik.grid <- apply(X.grid, 1, logLikFun, m)

# Plot optimal values from Theta 1 and Theta 2

contour(x.grid, x.grid, matrix(logLik.grid, n.grid, n.grid), 40, 
        xlab = expression(theta[1]),
        ylab = expression(theta[2]))

opt <- m@covariance@range.val

points(opt[1], opt[2], pch = 19, col = "red")

# Draw the Kriging Mean and vizualize accuracy

# Ploting the real function

n.grid <- 50

x.grid <- seq(0, 1, length = n.grid)

X.grid <- expand.grid(x1 = x.grid, x2 = x.grid)

y.grid <- apply(X.grid, 1, branin)

pred.m <- predict(m, X.grid, "UK")

par(mfrow = c(1, 3))

contour(x.grid, x.grid, matrix(y.grid, n.grid, n.grid), 50,
        main = "Braning")

# Plot the points to evaluete model

points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

# plot Kriging model

contour(x.grid, x.grid, matrix(pred.m$mean , n.grid, n.grid), 50,
        main = "Kriging Mean")

points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

# Plot Kriging Variance


contour(x.grid, x.grid, matrix(pred.m$sd^2 , n.grid, n.grid), 50,
        main = "Kriging Mean")

points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

plot(m)


# 4 Industrial Case with p = 6

rm(list = ls())

n <- 80 # number of points

d <- 6 # number of dimesions

set.seed(0)

X <- optimumLHS(n, d)

X <- data.frame(X)

y <- apply(X, 1, hartman6)

mlog <- km(design = X, response = -log(-y))

plot(mlog)

        plot.km(mlog)

n.test <- 250

set.seed(0)

X.test <- randomLHS(n.test, d)

colnames(X.test) <- names(X)

y.test <- apply(X.test, 1, hartman6)

ylog.pred <- predict(mlog, newdata = X.test, type = "UK")$mean

plot(ylog.pred)

# Sensitivity Analysis

branin.mat <- function(X) apply(X, 1, branin)

SA.branin <- fast99(model = branin.mat, 
                    factors = 2, 
                    n = 1000, 
                    q = "qunif",
                    q.arg = list(min = 0, max = 1))

par(mfrow = c(1, 2))

plot(SA.branin)

plot(SA.metamodel)

# Sensitivity with 8 dimension

n <- 80 # number of points

d <- 8 # number of dimesions

set.seed(0)

X <- optimumLHS(n, d)

X <- data.frame(X)

y <- sobol.fun(X)

m.sobol <- km(design = X, response = y)

SA.metamodel <- fast99(model = kriging.mean, 
                    factors = d, 
                    n = 1000, 
                    q = "qunif",
                    q.arg = list(min = 0, max = 1),
                    m = m.sobol)

SA.sobol <- fast99(model = sobol.fun, 
                    factors = d, 
                    n = 1000, 
                    q = "qunif",
                    q.arg = list(min = 0, max = 1))

par(mfrow = c(1, 2))

# Know problems

X <- expand.grid(x1 = seq(0, 1, length = 10), x2 = seq(0, 1, length = 10))
y <- branin(X)
t <- try(km(design = X, response = y, covtype = "gauss"))
cat(t)

t <- try(km(design = X, response = y, covtype = "gauss", nugget = 1e-8 * var(y)))

t.matern <- km(design = X, response = y)

