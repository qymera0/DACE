# 0 LOAD PACKAGES ---------------------------------------------------------

library(DiceEval)

library(DiceKriging)

library(DiceView)

library(tidyverse)

library(googlesheets)

library(leaps)

library(sensitivity)

library(VCA)

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


# 2 ORDINARY KRIGING ------------------------------------------------------

## 2.1 Wrangling Data

ok.exp <- transpose(cross(k.par))

## 2.2 Defining Covariance Structure with Cross Validation 

r.ok.sel <- ok.exp %>%
  pmap(modelFit, 
       X = X, 
       type = "Kriging", 
       multistart = 1,
       control = list(trace = FALSE)) 

r.ok.kfold <- r.ok.sel %>%
  map(crossValidation, K = dim(X)[1]) # LOO Rsquare

# Create dataframe with all LOO RSquares for all possibilities

r.ok.r2 <- setNames(data.frame(rep(names(k.par$Y), times = 10), 
                               rep(names(k.par$covtype), each = 4), 
                               rep(names(k.par$nugget.estim), each = 20),
                               stringsAsFactors = FALSE),
                    names(k.par))

r.ok.r2$R2 <- setNames(as.vector(lapply(r.ok.kfold, "[[", 2), mode = "numeric"),
                       c("R2"))
# Plot R Square Results

varPlot(form = as.formula(paste(names(r.ok.r2[4]), paste(names(r.ok.r2[-4]), collapse = " + "), sep = " ~ ")),
        Data = r.ok.r2,
        ylim = c(0.5, 1),
        Title = list("K Fold Performance for Ordinary Kriging"),
        YLabel = list(text = "LOO R Square"),
        MeanLine = list(var = "int"))

## 2.3 Final Ordinary Kriging Fit 

# sort and Filter the best configuration model

r.ok.best <- r.ok.r2 %>%
  group_by(Y) %>%
  filter(rank(desc(R2)) == 1) %>%
  arrange(Y)

# Create a list with all variables to allow model fit

k.par.final <- list(Y = list(r1 = y$r1, r2 = y$r2, r3 = y$r3, v7 = y$v7), 
                    covtype = list(r.ok.best$covtype[1],
                                   r.ok.best$covtype[2],
                                   r.ok.best$covtype[3],
                                   r.ok.best$covtype[4]),
                    nugget.estim = list(as.logical(r.ok.best$nugget.estim[1]),
                                        as.logical(r.ok.best$nugget.estim[2]),
                                        as.logical(r.ok.best$nugget.estim[3]),
                                        as.logical(r.ok.best$nugget.estim[4])))
# Fit Ordinary Kriging

r.ok <- k.par.final %>%
  pmap(modelFit, 
       X = X, 
       type = "Kriging",
       control = list(trace = FALSE),
       multistart = 10
  ) 

r.ok.kfold.final <- r.ok %>%
  map(crossValidation, K = dim(X)[1])
