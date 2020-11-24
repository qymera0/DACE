# 3.1 Discover Trend ------------------------------------------------------

x.trend <- X %>%
  select("Spin", "B_load", "U_load", "U_Height")

max.trend.models <- 14


# 3.1.1 K-Fold Forward Selection ------------------------------------------

x.trend2 <- x.trend

x.trend2$y1 <- y$r1

resp2 <- deparse(substitute(y1)) # Get the variable name as text

r1.pol.form2 <- as.formula(paste(resp2, 
                                 paste(".^2 + I(", paste(names(x.trend), collapse = "^2) + I("), "^2)"), 
                                 sep = " ~ ")) # Create the RSM formula for all variables


# Setting up K-Fold Cross Validation --------------------------------------

k <- 10 # Number of folds to make the selection

set.seed(1) # Set random seeds to replicable results

folds <- sample(1:k, nrow(x.trend2), replace = TRUE)

cv.errors <- matrix(NA, k, (max.trend.models - 4), 
                    dimnames = list(NULL, paste(1:(max.trend.models - 4)))) #K-Fold error matrix

r2.errors <- matrix(NA, k, (max.trend.models - 4), 
                    dimnames = list(NULL, paste(1:(max.trend.models - 4)))) #K-fold ^2 matrix

# Prediction method for forward selection

pred.kfwrd <- function (object, newdata, id, ...){
  
  form = as.formula(eval(object$call[[2]]))
  
  mat = model.matrix(form, newdata)
  
  coefi = coef(object, id = id)
  
  xvars = names(coefi)
  
  mat[ ,xvars]%*%coefi # Matrix multiplication
}

# K-Fold Forward Selection

for(j in 1:k) {
  
  r1.kffwrd.sel <- regsubsets(x = r1.pol.form2, 
                              data = x.trend2[folds != j, ],
                              method = "forward",
                              force.in = c(1, 2, 3, 4),
                              nvmax = max.trend.models 
  )
  
  for(i in 1:(max.trend.models - 4)){
    
    pred <- pred.kfwrd(r1.kffwrd.sel, x.trend2[folds == j, ], id = i)
    
    cv.errors[j, i] <- mean((x.trend2$y1[folds == j] - pred)^2)
    
    r2.errors[j, i] <- (cor(x.trend2$y1[folds == j], pred))^2
    
  } # max.trend.models - 3 because three elements are forced to enter
  
}

# K-Fold Model Plot

mean.cv.errors <- apply(cv.errors, 2, mean)

mean.r2 <- apply(r2.errors, 2, mean)

par(mfrow = c(2, 1))

plot(mean.cv.errors, 
     type = "b",
     main = "Cross Validation Errors",
     xlab = "Elements in the model",
     ylab = "Mean of Square Error")

points (which.min(mean.cv.errors), mean.cv.errors[which.min(mean.cv.errors)], col ="red", cex =2, pch =20)

plot(mean.r2, 
     type = "b",
     main = "Cross Validation R Square",
     xlab = "Elements in the model",
     ylab = "K Fold R Square")

points (which.max(mean.r2), mean.r2[which.max(mean.r2)], col ="red", cex =2, pch =20)

# Final trend Model Adjust

resp <- deparse(substitute(y$r1)) # Get the variable name as text

r1.pol.form <- as.formula(paste(resp, 
                                paste(".^2 + I(", paste(names(x.trend), collapse = "^2) + I("), "^2)"), 
                                sep = " ~ ")) # Create the RSM formula for noise variables

r1.trend <- regsubsets(x = r1.pol.form, 
                       data = x.trend,
                       method = "forward",
                       force.in = c(1, 2, 3, 4))                                           


# Predict Using Final trend model

r1.trend.form <- as.formula(paste(resp, 
                                  paste(paste0(names(coef(r1.trend, 
                                                          which.min(mean.cv.errors))[2:(which.min(mean.cv.errors) + 5)])), 
                                        collapse = " + "),
                                  sep = " ~ "))


r1.trend.m <- lm(r1.trend.form)

summary(r1.trend.m)

par(mfrow = c(1,1))

plot(predict(r1.trend.m), y$r1,
     main = "Predictions for 10 Fold UK Trend",
     ylab = "Radial Displacement at Balance Ring (rad1)",
     xlab = "Trend Prediction for rad1")


# 3.2 Fit Universal Kriging -----------------------------------------------

r1.uk <- km(formula = r1.trend.form,
            design = X, 
            response = y$r1,
            covtype = "gauss",
            nugget.estim = TRUE,
            optim.method = "gen")


plot.km(r1.uk)


# 3.3 Sensitivity Analysis ------------------------------------------------

ranges <- list()

for (i in names(X)) {
  
  ranges[[i]] <- list(min(eval(parse(text = paste("X$", i)))), 
                      max(eval(parse(text = paste("X$", i)))))
  
  names(ranges[[i]]) <- c("min", "max")
  
}

r1.mean <- function(Xnew, m) predict.km(m, Xnew, "UK", checkNames = "FALSE")$mean

r1.sens <- fast99(model = r1.mean, 
                  factors = 11,
                  n = 1000,
                  q = "qunif",
                  q.arg = ranges,
                  m = r1.OK)

names(r1.sens$X) <- names(X)

par(mfrow = c(1,1))

plot(r1.sens, 
     main = "Rad1 Ordinary Kriging Variable Importance",
     xlab = "Variable",
     ylab = "Standardized Importance")

