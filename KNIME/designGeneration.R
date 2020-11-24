library(DiceDesign)

# 01 INITIAL CONFIGURATION ------------------------------------------------

m <- knime.flow.in[["iterations"]] # number of optimization repetitions

n <- knime.flow.in[["nruns"]] # number of runs

p <- length(knime.in$"Factor") # number of factors

factors <- knime.in

# 02 INITIAL LHS DESIGN ---------------------------------------------------

lhs <- lhsDesign(n, p, seed = 123456)

# 03 DESIGN REPLICATES ----------------------------------------------------

# Initiate an empty list

design <- list()

# Half of designs uses genetic algorithm

for(i in 1:(m %/%2)){
        
        design[[i]] <- discrepESE_LHS(lhs$design, 
                                      criterion = "L2star",
                                      inner_it = 100,
                                      J = 50,
                                      it = 5)
        
}

# Half of designs uses simulated annealing 

for(i in ((m %/%2)+1):m){
        
        design[[i]] <- discrepSA_LHS(lhs$design)
        
}

# 04 CALCULATE RSS FOR ALL DESIGNS ----------------------------------------

# Initialize rss list

rssAll <- list()

# Calculate all rss statistics for all design

for(i in 1:m){
        
        rssAll[[i]] <- rss2d(design[[i]]$design,
                             lower = seq(from = 0,
                                         to = 0,
                                         length.out = p),
                             upper = seq(from = 1,
                                         to = 1,
                                         length.out = p),
                             graphics = -1,
                             trace = F)
        
}

# 05 SELECT A DESIGN ------------------------------------------------------

# Create a vector with the differences between max and limit RSS

# Initialize vector

maxRss <- vector()

# Create a dataframe with the differences between max and limit RSS

for(i in 1:m) {
        
        maxRss[i] <- rssAll[[i]]$gof.test.stat - max(rssAll[[i]]$global.stat, 
                                                     na.rm = T)      
        
}

# 06 EXPORT DATA BACK ------------------------------------------------------

finallDesign <- design[[which.max(maxRss)]] # list with all design infomration

knime.out <- as.data.frame(unscaleDesign(finallDesign$design, 
                                         min = factors$"level.minus",
                                         max = factors$"level.plus")$design)

names(knime.out) <- knime.in$"Factor"

# 07 DESIGN QUALITY METRICS -----------------------------------------------

rss <- rss2d(knime.out,
             lower = factors$level.minus,
             upper = factors$level.plus,
             graphics = -1,
             trace = F)

knime.out <- data.frame(Metric = c("coverage", 
                                   "meshRatio", 
                                   "mindist",
                                   "discC2",
                                   "discL2",
                                   "discM2",
                                   "discS2",
                                   "discW2",
                                   "rssDiff"),
                        Value = c(coverage(knime.out),
                                  meshRatio(knime.out)[[1]],
                                  mindist(knime.out),
                                  discrepancyCriteria(knime.out, type = "C2")[[1]],
                                  discrepancyCriteria(knime.out, type = "L2")[[1]],
                                  discrepancyCriteria(knime.out, type = "M2")[[1]],
                                  discrepancyCriteria(knime.out, type = "S2")[[1]],
                                  discrepancyCriteria(knime.out, type = "W2")[[1]],
                                  rss$gof.test.stat - max(rss$global.stat, na.rm = T)))