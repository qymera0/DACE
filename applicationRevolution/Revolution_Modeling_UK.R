# 0 LOAD PACKAGES ---------------------------------------------------------

library(DiceEval)

library(DiceKriging)

library(DiceView)

library(DiceDesign)

library(tidyverse)

library(googlesheets)

library(leaps)

library(sensitivity)

library(VCA)

library(triangle)

library(pse)

# 1 LOAD DATA -------------------------------------------------------------

full.dace.key <- gs_key("1pGFVdeUy2kqAiy8KoCiANibrV0rrYN7PqrA-lHJvHt0", 
                        visibility = "private") # Get Sheet for GDrive

brz <- locale(decimal_mark = ",") # Define "," as decimal mark

full.dace <- gs_read(full.dace.key, ws = "Big", locale = brz)

X <- full.dace %>%
  select(-"r1", -"r2", -"r3", -"v7")

y <- full.dace %>%
  select("r1", "r2", "r3", "v7")

rm(brz, full.dace.key)

set.seed(123456)

# Create a list with all Kriging fit parameters

k.par <- list(Y = list(r1 = y$r1, r2 = y$r2, r3 = y$r3, v7 = y$v7), 
              covtype = list(gauss = "gauss", 
                             matern5_2 = "matern5_2", 
                             matern3_2 = "matern3_2", 
                             exp = "exp", 
                             powexp = "powexp"), 
              nugget.estim = list(True = TRUE, False = FALSE))


# 2 WRANGLING DATA --------------------------------------------------------

x.trend <- X %>%
  select("Spin", "B_load", "U_load", "U_Height")

x.factors <- X %>%
  select(-"Spin", -"B_load", -"U_load", -"U_Height")

x.factors$Spin <- X$Spin


rsm.trend.for <- paste0(".^2 + I(", paste(names(x.trend), collapse = "^2) + I("), "^2)")

uk.trend.par <- list(Y = list(r1 = y$r1, r2 = y$r2, r3 = y$r3, v7 = y$v7),
                     formula = list(r1 = paste("Y", rsm.trend.for, sep = " ~ "),
                                    r2 = paste("Y", rsm.trend.for, sep = " ~ "),
                                    r3 = paste("Y", rsm.trend.for, sep = " ~ "),
                                    v7 = paste("Y", rsm.trend.for, sep = " ~ ")
               )
)



# 3 DISCOVER TREND --------------------------------------------------------

# Define the best penalty constant for stepwise trend formula

uk.trend.sel <- uk.trend.par %>%
  pmap(stepEvolution,
       X = x.trend,
       P = 1:30,
       graphic = FALSE
       ) 

# Complement the list of parameters with optimal penalty parameters

uk.trend.par$penalty <- list(penalty = which.max(uk.trend.sel$r1$Q2),
                             penalty = which.max(uk.trend.sel$r2$Q2),
                             penalty = which.max(uk.trend.sel$r3$Q2),
                             penalty = which.max(uk.trend.sel$v7$Q2))

# Make optimal forward selection for trend formula

uk.trend <- uk.trend.par %>%
  pmap(modelFit,
       X = x.trend,
       type = "StepLinear")

# 4 DISCOVER OPTIMAL MODEL CONFIGURATION ----------------------------------

# Wrangling data to UK configurations experiment

uk.exp <- transpose(cross(k.par))

uk.exp$formula <- rep(list(uk.trend$r1$model$terms,
                           uk.trend$r2$model$terms,
                           uk.trend$r3$model$terms,
                           uk.trend$v7$model$terms), times = 10)

# Fit UK Model for all possible configurations

r.uk.sel <- uk.exp %>%
  pmap(modelFit, 
       X = X, 
       type = "Kriging", 
       multistart = 1,
       control = list(trace = FALSE)) 

# K Fold for all possibilities using Leave One Out (LOO)

safe_cross <- possibly(crossValidation, list(0, 0))

r.uk.kfold <- r.uk.sel %>%
  map(safe_cross, K = dim(X)[1]) # LOO Rsquare

# Create dataframe with all LOO RSquares for all possibilities

r.uk.r2 <- setNames(data.frame(rep(names(k.par$Y), times = 10), 
                               rep(names(k.par$covtype), each = 4), 
                               rep(names(k.par$nugget.estim), each = 20),
                               stringsAsFactors = FALSE),
                    names(k.par))

r.uk.r2$R2 <- setNames(as.vector(lapply(r.uk.kfold, "[[", 2), mode = "numeric"),
                       c("R2"))

# Plot R Square Results

varPlot(form = as.formula(paste(names(r.uk.r2[4]), paste(names(r.uk.r2[-4]), collapse = " + "), sep = " ~ ")),
        Data = r.uk.r2,
        ylim = c(0.8, 1),
        Title = list("K Fold Performance for Universal Kriging"),
        YLabel = list(text = "LOO R Square"),
        MeanLine = list(var = "int"))

# Select the best configuration for each response

r.uk.best <- r.uk.r2 %>%
  group_by(Y) %>%
  filter(rank(desc(R2)) == 1) %>%
  arrange(Y)


# 5 FIT FINAL MODEL -------------------------------------------------------

# Create a list with all variables to allow model fit

k.par.final <- list(Y = list(y$r1, y$r2, y$r3, y$v7), 
                    covtype = list(r.uk.best$covtype[1],
                                   r.uk.best$covtype[2],
                                   r.uk.best$covtype[3],
                                   r.uk.best$covtype[4]),
                    nugget.estim = list(as.logical(r.uk.best$nugget.estim[1]),
                                        as.logical(r.uk.best$nugget.estim[2]),
                                        as.logical(r.uk.best$nugget.estim[3]),
                                        as.logical(r.uk.best$nugget.estim[4])),
                    formula = list(uk.trend$r1$model$terms,
                                   uk.trend$r2$model$terms,
                                   uk.trend$r3$model$terms,
                                   uk.trend$v7$model$terms))

# Fit Final Universal Kriging

r.uk <- k.par.final %>%
  pmap(modelFit, 
       X = X, 
       type = "Kriging",
       control = list(trace = FALSE),
       multistart = 10
  ) 

names(r.uk) <- names(y)

# Evaluate final models with K-Fold Cross Validation

r.uk.kfold.final <- r.uk %>%
  map(crossValidation, K = dim(X)[1])

# Create a model list

model.list <- list(m = r.uk$r1$m,
                   m = r.uk$r2$m,
                   m = r.uk$r3$m,
                   m = r.uk$v7$m)

# 6 QUALITY OF FIT --------------------------------------------------------

model.list %>%
  map(plot,
      trend.reestim = TRUE)

# 7 VARIABLES IMPORTANCE --------------------------------------------------

# Create a list with minimum and maximum for all variables

x.range <- vector("list", dim(X)[2])

for (i in 1:dim(X)[2]) {
  
  x.range[[i]] <- list(min = min(X[i]), max = max(X[i]))
  
}

names(x.range) <- names(X)

# Create other parameter for graphics

x.dist <- as.vector(rep("qunif", dim(X)[2]))

x.center <- colMeans(X)

# Plot Marginal Model Graphics

model.list %>%
  map(sectionview.km,
      center = x.center,
      conf_lev = 0.95,
      conf_blend = 0.2,
      bg_blend = 1,
      title = " ",
      Xname = names(x.range),
      yname = "comb. Rad. Disp [mm]",
  )

#Function to make predictions using kriging models

kriging.mean <- function(Xnew, m) predict.km(m, Xnew, "UK", se.compute = FALSE, checkNames = FALSE)$mean

# Sensitivity Analysis for Variables

par(mfrow = c(1,1))

sa.uk <- vector("list", dim(y)[2])

j <- 1

for (i in model.list) {
  
  sa.uk[[j]] <- fast99(model = kriging.mean,
                       factors = names(X),
                       n = 1000,
                       q = x.dist,
                       q.arg = x.range,
                       m = i)
  plot(sa.uk[[j]],
       i = 1:7)
  
  title(main = paste("Sensitivity Analysis for ", names(y[j])))
  
  j = j + 1
  
  
}

rm(i, j)

names(sa.uk) <- names(y)

# 8 SIMULATE TO modeFRONTIER ----------------------------------------------

mc.size <- 10000

exp.runs <- 100

# modeFRONTIER FIT

# Parameters to simulate noise factors

noise.par <- list(B_load = list(a = 0, b = 22, c = 4.94975),
                  U_load = list(a = 0, b = 2.5, c = 0.70711),
                  U_Height = list(a = 50, b = 275, c = 141.421))


# Simulate Noise Factors

noise.sim <- transpose(noise.par) %>%
  pmap(rtriangle, n = mc.size) %>%
  as.data.frame() %>%
  slice(rep(1:n(), times = exp.runs))

# Create a list with minimum and maximum for design variables

x.des.r <- vector("list", dim(x.factors)[2])

for (i in 1:dim(x.factors)[2]) {
  
  x.des.r[[i]] <- list(min = min(x.factors[i]), max = max(x.factors[i]))
  
}

names(x.des.r) <- names(x.factors)

# Create a list with factor distributions

x.des.dist <- as.vector(rep("qunif", dim(x.factors)[2]))

# Create the experimental desing for design variables

x.lhs <- LHS(factors = names(x.factors),
             N = exp.runs,
             q = x.des.dist,
             q.arg = x.des.r)$data %>%
  slice(rep(1:n(), each = mc.size))

x.sim <- cbind(x.lhs, noise.sim)

d <- dim(x.sim)[2]

for (i in 1:dim(y)[2]) {
  
  pi <- paste0('p', names(y[i]))
  
  x.sim[ ,pi] <- predict.km(model.list[i], x.sim[1:d], "UK", se.compute = FALSE, checkNames = FALSE)
  
  
}
