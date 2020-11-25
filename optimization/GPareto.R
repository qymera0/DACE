library(DiceDesign)
library(DiceEval)
library(DiceKriging)
library(DiceOptim)
library(DiceView)
library(GPareto)

# TWO OBJECITVES, UNIDIMENSIONAL ------------------------------------------

# Model fitting

design.init <- matrix(seq(0, 1, length.out = 6), ncol = 1) # Initial Design

response.init <- MOP2(design.init) # Initial response evaluation

mf1 <- km(~1, design = design.init, response = response.init[, 1])

mf2 <- km(~1, design = design.init, response = response.init[, 2])

model <- list(mf1, mf2)

# Optimization

res <-
        GParetoptim(
                model = model,
                fn = MOP2,
                crit = "EHI",
                nsteps = 7,
                lower = 0,
                upper = 1,
                critcontrol = list(refPoint = c(2,2))
        )

# MORE COMPLEX OPTIMIZATION -----------------------------------------------

# Model fitting

set.seed(1)

d <- 2; ninit <- 7; fun <- P1

design <- lhsDesign(ninit, d, seed = 42)$design

response <- t(apply(design, 1, fun))

mf1 <- km(~., design = design, response = response[, 1])

mf2 <- km(~., design = design, response = response[, 2])

model <- list(mf1, mf2)

# Optimization

x.grid <- seq(0, 1, length.out = 21)

test.grid <- expand.grid(x.grid, x.grid)

SURcontrol <- list(integration.points = test.grid)

omEGO1<-
        crit_optimizer(
                crit = "SUR",
                model = model,
                lower = c(0,0),
                upper = c(1,1),
                critcontrol = list(SURcontrol = SURcontrol),
                optimcontrol = list(method = "genoud", 
                                    pop.size = 20,
                                    int.seed = 2,
                                    unif.seed = 3)
        )

# Separate functions

fun1 <- function(x) P1(x)[ ,1]

fun2 <- function(x) P1(x)[ ,2]

fastmf2 <- fastfun(fn = fun2,
                   design = design,
                   response = response[ ,2])

model2 <- list(mf1, fastmf2)

# Optimization - 

sol <-
        GParetoptim(
                model = model,
                fn = fun,
                crit = "SUR",
                nsteps = 7,
                lower = c(0,0),
                upper = c(1,1),
                optimcontrol = list(method = "pso"),
                critcontrol = list(SURcontrol = list(distrib = "SUR", 
                                                     n.points = 50)
                )
        ) # demora pra caramba

solFast <-
        GParetoptim(
                model = list(mf1),
                fn = fun1,
                cheapfn = fun2,
                crit = "SUR",
                nsteps = 7,
                lower = c(0,0),
                upper = c(1,1),
                optimcontrol = list(method = "pso"),
                critcontrol = list(SURcontrol = list(distrib = "SUR", 
                                                     n.points = 50)
                )
        )
