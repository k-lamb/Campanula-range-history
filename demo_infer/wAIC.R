# build boxplot of model outcomes
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Documents/Research/Paper Code/Debban-Lamb/FINAL/demo_infer/data_concat/")

perc <- 65 # change to percentage retained to call correct files (80/65)
meta <- F # working on population aggregates (TRUE) or individual populations (FALSE)

# read in data and add columns where necessary to allow data frames to merge correctly
if (meta == T) {
  if (perc == 80) {
    dat1 <- read.csv("AIC_concat_80perc_linreg.csv") # everything. change for differences in filtering (656perc vs 80perc) and lineage and region compositing (linreg)
    dat2 <- read.csv("AIC_concat_80perc_linreg_2c.csv")
  }
  if (perc == 65) {
    dat1 <- read.csv("AIC_concat_65perc_linreg.csv") # everything. change for differences in filtering (656perc vs 80perc) and lineage and region compositing (linreg)
    dat2 <- read.csv("AIC_concat_65perc_linreg_2c.csv")
  }
  # create empty/NA variables so dat1 and dat2 can be merged
  dat1$Ts2 <- NA
  dat1$Tsc2 <- NA
  dat1$m12_2 <- NA
  dat1$m21_2 <- NA
  dat1 <- dat1[, names(dat2)] # get columns in same order in dat1 as in dat2
  dat <- rbind(dat1,dat2)
  test <- dat
}

if (meta == F) {
  if (perc == 80) {
    dat1 <- read.csv("AIC_concat_80perc.csv") # singular data
    dat2 <- read.csv("AIC_concat_80perc_2c.csv") # cyclic data
    
    # create empty/NA variables so dat1 and dat2 can be merged
    dat1$Ts2 <- NA
    dat1$Tsc2 <- NA
    dat1$m12_2 <- NA
    dat1$m21_2 <- NA
    dat1 <- dat1[, names(dat2)] # get columns in same order in dat1 as in dat2
    dat <- rbind(dat1,dat2)
    test <- dat
  }
  if (perc == 65) {
    test <- read.csv("AIC_concat_65perc.csv")
  }
}

# dat1 <- read.csv("AIC_concat_80perc_linreg.csv") # everything. change for differences in filtering (656perc vs 80perc) and lineage and region compositing (linreg)
# dat2 <- read.csv("AIC_concat_80perc_linreg_2c.csv")


# to see only the main models
# test <- test %>% filter(model == "IM_ae" | model == "SC_ae" | model == "NM_ae" | model == "AM_ae")
# test <- test %>% filter(model == "IM" | model == "SC" | model == "NM" | model == "AM" | model == "A_SC" |
#                         model == "A_SC" | model == "IM2C" | model == "SC2C")

# calculate akaike weights
aic.min <- aggregate(aic~Pair_name, data=test, FUN="min") # find min aic for pop pair
names(aic.min)[2] <- "aic.min" # min score for delta(aic) calculation
t <- merge(test, aic.min, all=T) 
t$aic.delta <- t$aic - t$aic.min

# t <- subset(t, aic.delta < 10) # remove highly improbable models

t$aic.d.exp <- exp((-t$aic.delta)/2) # numerator of aic.w calculation as per wagenmakers and farrell 2004
aic.sum <- aggregate(aic.d.exp~Pair_name, data=t, FUN="sum") # denominator of aic.w equation
names(aic.sum)[2] <- "aic.d.exp.sum"
t <- merge(t, aic.sum, all=T)

t$aic.w <- t$aic.d.exp/t$aic.d.exp.sum # weighted aic

# calculate aic weights with MuMIn to see comparison
# matches my script exactly :)
# t <- t %>% 
#   group_by(Pair_name) %>% 
#   mutate(aic.mumin = MuMIn::Weights(aic))

# ratio of weights to see likelihood -- ignore
aic.w.best <- aggregate(aic.w~Pair_name, data=t, FUN="max") # best weight per pair
names(aic.w.best)[2] <- "aic.w.best"
t <- merge(t, aic.w.best, all=T)
t$aic.ll <- t$aic.w.best/t$aic.w

# adding metric of model quality per Rougeux et al. 2017

# find maximum ∆AIC
aic.d.max <- aggregate(aic.delta~Pair_name, data=t, FUN="max")
names(aic.d.max)[2] <- "aic.d.max"
t <- merge(t, aic.d.max, all=T)

# model quality index calculation
t$model_score <- (t$aic.d.max - t$aic.delta) / t$aic.d.max

# set model score NaN to 1 (I think this comes from being uncontested, so AIC d max = 0)
t <- t %>% mutate(model_score = case_when(
  model_score != "NaN" ~ model_score,
  model_score == "NaN" ~ 1
  )
)

# write to csv
if (meta == TRUE) {
  if (perc == 80) {
    write.csv(t, "AIC_weights_all_80perc_linreg_ext_2c_ALL.csv", row.names = F)
  }
  if (perc == 65) {
    write.csv(t, "AIC_weights_all_65perc_linreg_ext_2c_ALL.csv", row.names = F)
  }
}

if (meta == FALSE) {
  if (perc == 80) {
    write.csv(t, "AIC_weights_all_80perc.csv", row.names = F)
  }
  if (perc == 65) {
    write.csv(t, "AIC_weights_all_65perc.csv", row.names = F)
  }
}

