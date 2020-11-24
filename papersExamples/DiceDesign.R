# Dice Design

library(DiceDesign)

n <- 20 # Number of runs

dimension <- 2 # number of variables

X_random <- matrix(runif(n * dimension), ncol = dimension, nrow = n) # generate a random design

x <- seq(0, 1, length = 5) # create a sequence of 5 equaly distributed numbers

grid <- expand.grid(x, x) # create a grid of equally spaced values

coverage(X_random)

meshRatio(X_random)

mindist(X_random)

discrepancyCriteria(X_random, type = c("L2", "M2"))

rss <- rss2d(design = X_random, lower = rep(0, dimension), upper = rep(1, dimension))

rss2 <- rss2d(design = grid, lower = rep(0, dimension), upper = rep(1, dimension))

# Spacefilling designs

X1 <- dmaxDesign(n = 20, dimension = 2, range = 0.2, niter_max = 1000, seed = 1)

dim <- 2

n <- 20

X2 <- straussDesign(n, dim, RND = 0.21, alpha = 0.2, repulsion = 0.0001)

# Uniform Optimal Space Filling Designs

n <- 30

dim <- 12

X <- lhsDesign(n, dim)$design # Initial Design

x_unif <- discrepESE_LHS(X, criterion = "L2star")

coverage(x_unif$design)

meshRatio(x_unif$design)

mindist(x_unif$design)

discrepancyCriteria(x_unif$design, type = "L2star")

rss <- rss2d(design = x_unif$design, lower = rep(0, dim), upper = rep(1, dim))
