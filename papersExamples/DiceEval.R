library(DiceEval)

# Fitting a model

data("dataIRSN5D")

pairs(dataIRSN5D, bg = "black", pch = 19)

pairs(testIRSN5D, bg = "black", pch = 19)

X <- dataIRSN5D[ ,-6]

Y <- dataIRSN5D[ ,6]

data("testIRSN5D")

X.test <- testIRSN5D[ ,-6]

Y.test <- testIRSN5D[ ,6]

test.data <- X.test

test.data$Y <- Y.test

modeLm <- modelFit(X, Y, type = "Linear", formula = Y ~ .)

names (modeLm)

summary(modeLm$model)

kfold <- crossValidation(modeLm, K=10)

residualsStudy(modeLm)

library("gam")

modeAm <- modelFit(X, Y, type = "Additive", formula = formulaAm(X, Y))

residualsStudy(modeAm)

modeStep <- modelFit(X, Y, type = "StepLinear",
                     formula = Y ~ .^2 + I(e^2) + I(r^2) + I(p^2) + I(l^2) + I(b^2),
                     penalty = 2)

modeStep2 <- modelFit(X, Y, type = "Linear",
                     formula = Y ~ .^2,
                     penalty = 2)

library("mda")

modMARS <- modelFit(X, Y, type = "MARS", degree = 2)

library("polspline")

modPolyMARS <- modelFit(X, Y, type = "PolyMARS", gcv = 4)

# Evaluating the Penalty Parameter at Stepwise

out <- stepEvolution(X, Y,
                     Y ~ .^2 + I(e^2) + I(r^2) + I(p^2) + I(l^2) + I(b^2),
                     P = c(1, 2, 3, 4, 5, 10, 20, 30))



out$Q2

# Evaluating the penalty for Polymars

Crit <- penaltyPolyMARS(X, Y, graphic = TRUE)

# model comparison

crit <- modelComparison(X, Y,
type = c("Linear", "Additive", "MARS", "PolyMARS", "Kriging"),
test = test.data,
degree = 2,
gcv = 2,
formula = c(Y ~ .^2 + I(e^2) + I(r^2) + I(p^2) + I(l^2) + I(b^2),
formulaAm(X, Y), Y ~ .))


