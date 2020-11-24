# 0 LOAD PACKAGES ---------------------------------------------------------

library(DiceEval)
library(DiceView)
library(tidyverse)
library(readr)
library(furrr)
library(magrittr)
library(doParallel)
library(sensitivity)

# set multi process environment

plan(multiprocess)

set.seed(123456)

# 1 LOAD DATA -------------------------------------------------------------

X <-
        knime.in %>%
        select(str_squish(unlist(str_split(knime.flow.in[["Concatenate (#1)"]], ","))))

y <-
        knime.in %>%
        select(str_squish(unlist(str_split(knime.flow.in[["Concatenate"]], ","))))

# 2 KRIGIN PARAMETERS -----------------------------------------------------

kPar <- list(Y = t(as.list(y)), 
             covtype = list(gauss = "gauss", 
                            matern5_2 = "matern5_2", 
                            matern3_2 = "matern3_2", 
                            exp = "exp", 
                            powexp = "powexp"), 
             nugget.estim = list(True = TRUE, False = FALSE))

# 3 TREND DETERMINATION ---------------------------------------------------

# Trend formula

rsmTrendFor <- paste0(".^2 + I(", paste(names(X), collapse = "^2) + I("), "^2)")

# Trend parameters

trendPar <- list(Y = kPar[[1]])

for(i in 1:length(y)) trendPar$formula[[i]] <- paste("Y", rsmTrendFor, sep = " ~ ")

# Trend identification

trendSel <- 
        trendPar %>%
        future_pmap(stepEvolution,
                    X = X,
                    P = 1:30,
                    graphic = FALSE)

# Adding optimal penalty parameters

for(i in 1:length(y)) trendPar$penalty[[i]] <- which.max(trendSel[[i]]$Q2)

# Make the final trend fit

trend <-
        trendPar %>%
        future_pmap(modelFit,
                    X = X,
                    type = "StepLinear")

# 4 KRIGIN MODEL FIT ------------------------------------------------------

# Prepare list of possibilities

kParDoe <- transpose(cross(kPar)) # All possible configurations

# Add trend formulas

kParDoe$formula <- rep((trend %>%
                                map(extract2, "model") %>%
                                map(extract2, "terms")),
                       times = 10)

# Fit Universal Kriging to all possibilities

kSel <-
        kParDoe %>%
        future_pmap(modelFit, 
                    X = X, 
                    type = "Kriging", 
                    multistart = 1,
                    control = list(trace = FALSE))

# Select the configurations that have the bigger Q2 (kfold R2)

safeCross <- possibly(crossValidation, list(0, 0)) # Capture possible errors

kKfoldMaxQ2 <- 
        kSel %>%
        future_map(safeCross, 
                   K = dim(X)[1]) %>%
        map(magrittr::extract, "Q2") %>%
        bind_rows() %>%
        mutate(idx = seq(from = 1, to = length(y)*10)) %>%
        mutate(y = rep(seq(from = 1, 
                           to = length(y)), times = 10)) %>%
        group_by(y) %>%
        filter(Q2 == max(Q2)) %>%
        distinct(Q2, .keep_all = T)

# 5 FINAL MODEL FIT -------------------------------------------------------

# Create a list with best configurations

kParFinal <- list(Y = kPar[[1]])

for(i in kKfoldMaxQ2$idx){
        kParFinal$covtype[i] <- kParDoe$covtype[i]
        kParFinal$nugget.estim[i] <- kParDoe$nugget.estim[i]
        kParFinal$formula[i] <- kParDoe$formula[i]
} 

# Clean NULL elements

kParFinal <-
        kParFinal %>%
        map(discard, is.null)

names(kParFinal$Y) <- paste(kKfoldMaxQ2$y)
names(kParFinal$covtype) <- paste(kKfoldMaxQ2$y)
names(kParFinal$nugget.estim) <- paste(kKfoldMaxQ2$y)
names(kParFinal$formula) <- paste(kKfoldMaxQ2$y)

# Fit final Universal Kriging

workers <- detectCores() - 1

cl <- makeCluster(workers)

registerDoParallel(cl)

uk <-
        kParFinal %>%
        pmap(modelFit, 
             X = X, 
             type = "Kriging",
             control = list(trace = FALSE),
             multistart = 30)

stopCluster(cl)

registerDoSEQ()

uk <- uk[order(names(uk))]

names(uk) <- names(y)

modelList <- list()

for(i in 1:length(y)) {
        
        modelList[[i]] <- uk[[i]]$model
        
}

names(modelList) <- names(y)


# PLOT VARIABLE IMPORTANCE ------------------------------------------------


k_mean <- function(Xnew, m) predict.km(m, 
                                       Xnew, 
                                       "UK", 
                                       se.compute = FALSE, 
                                       checkNames = FALSE)$mean
par(mfrow = c(1,1))

# List with variables range

sensList <- list(xDist = as.vector(rep("qunif", dim(X)[2])),
                 xRange = vector("list", dim(X)[2])) 

for (i in 1:dim(X)[2])  sensList$xRange[[i]] <- list(min = min(X[i]), 
                                                     max = max(X[i]))

names(sensList$xRange) <- names(X)

# Sensitivity Analysis

saUk <- list()

j <- 1

for (i in modelList) {
        
        saUk[[j]] <- fast99(model = k_mean,
                            factors = names(X),
                            n = 1000,
                            q = sensList$xDist,
                            q.arg = sensList$xRange,
                            m = i)
        
        plot(saUk[[j]], xprob = T)
        
        title(main = paste("Sensitivity Analysis for ", names(y)[j]))
        
        png(paste(knime.flow.in[["context.workflow.absolute-path"]],"\\sens", names(modelList[j]),".png", sep = ""))
        
        plot(saUk[[j]], xprob = T)
        
        title(main = paste("Sensitivity Analysis for ", names(y)[j]))
        
        dev.off()
        
        j = j + 1
        
}

# PLOT RESIDUALS ----------------------------------------------------------

for (i in 1:length(modelList)){
        
        plot(modelList[[i]])
        title(names(modelList)[i], adj = 0)
        png(paste(knime.flow.in[["context.workflow.absolute-path"]],"\\residual", names(modelList[i]),".png", sep = ""))
        plot(modelList[[i]])
        dev.off()
        
}


# PLOT SECTION VIEW -------------------------------------------------------

# Plot Section View
for (i in 1:length(modelList)){
        
        sectionview(modelList[[i]],
                    center = colMeans(X),
                    conf_lev = 0.95,
                    conf_blend = 0.2,
                    bg_blend = 1,
                    Xname = names(X),
                    yname = names(modelList[i]),
                    title = " ")
        png(paste(knime.flow.in[["context.workflow.absolute-path"]],"\\section", names(modelList[i]),".png", sep = ""))
        sectionview(modelList[[i]],
                    center = colMeans(X),
                    conf_lev = 0.95,
                    conf_blend = 0.2,
                    bg_blend = 1,
                    Xname = names(X),
                    yname = names(modelList[i]),
                    title = " ")
        dev.off()
        
}



