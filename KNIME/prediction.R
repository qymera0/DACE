# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(magrittr)
library(DiceDesign)
library(DiceEval)
library(DiceKriging)
library(DiceView)

# 1 LOAD DATA -------------------------------------------------------------

ccd <- knime.in

# 2 LOAD R MODEL ----------------------------------------------------------

ccdPredict <-
        modelList %>% 
        map(predict, 
            newdata = ccd,
            type = "UK",
            light.return = TRUE,
            se.compute = FALSE) %>%
        map(extract2, "mean") %>%
        as_tibble()

knime.out <- cbind(ccd, ccdPredict)


