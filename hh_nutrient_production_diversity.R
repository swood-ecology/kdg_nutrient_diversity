library(haven)    # For working with Stata data
library(RStata)   # For working with Stata data
library(FD)       # For calculating FD metrics
source("~/Box Sync/Work/Statistics/Code/setRows.R")  # for setting rows

# Production data
RICARDIANDATA <- read_dta("~/Box Sync/Work/Writing/Manuscripts/Submitted/Fonio nutritional diversity/RICARDIANDATA.dta")

# Nutrient data
nutr <- read.csv("~/Box Sync/Work/Writing/Manuscripts/Submitted/Fonio nutritional diversity/nutrientData.csv")

# Create new data frame of production data in kg
rice <- RICARDIANDATA$riceSacks*50
peanuts <- RICARDIANDATA$peanutsSacs*50
millet <- RICARDIANDATA$milletSacks*50
maize <- RICARDIANDATA$cornSacks*50
fonio <- RICARDIANDATA$fonioSacs*50
taro <- RICARDIANDATA$tarotSacks*50
potatoes <- RICARDIANDATA$potatoAmtHarvested

hh <- data.frame(rice,peanuts,millet,maize,fonio,taro,potatoes)
hh[hh == 0] <- NA
hh[126,4] <- NA

# Convert nutrient data to be per kg
nutr <- set.rownames(nutr)
rownames(nutr) <- c('fonio','maize','rice','peanuts','sorghum','cowpea','cassava','cabbage','potatoes','okra','bissap','jakatu','millet','taro')
nutr <- nutr*10
nutr <- nutr[,-c(13,15,22)]

# Crop species with nutritient data but not found in farms. If blank, then all is well.
rownames(nutr)[!rownames(nutr) %in% colnames(hh)]

# Drop unobserved crops from nutrient database
nutr <- nutr[ ! rownames(nutr) %in% rownames(nutr)[!rownames(nutr) %in% colnames(hh)], ]

# Order data frames for FD metrics
nutr <- nutr[order(rownames(nutr)),]
hh <- hh[,order(colnames(hh))]

# Calculate diversity
diversity_scores <- dbFD(nutr,hh,calc.CWM=F,calc.FGR=T) # 4 groups
diversity_data_frame <- data.frame(as.factor(diversity_scores$nbsp),diversity_scores$FDis,diversity_scores$FDiv,diversity_scores$FEve,diversity_scores$FRic)
names(diversity_data_frame) <- c("SpeciesNumber","FDis","FDiv","FEve","FRic")
pairs(diversity_data_frame[,-1], upper.panel=NULL)

par(mfrow=c(2,2))
boxplot(FDis~SpeciesNumber,data=diversity_data_frame,ylab="Functional dispersion (FDis)")
boxplot(FDiv~SpeciesNumber,data=diversity_data_frame,ylab="Functional divergence (FDiv)")
boxplot(FEve~SpeciesNumber,data=diversity_data_frame,ylab="Functional evenness (FEve)")
boxplot(FRic~SpeciesNumber,data=diversity_data_frame,ylab="Functional richness (FRic)")

# Simulations
diversity_simulations <- simul.dbFD(s=c(3, 5, 10, 15),t=20, r=200, p=15, 
                                    tr.method="unif", abun.method="lnorm", 
                                    w.abun=T)
diversity_simulations_low_traits <- simul.dbFD(s=c(3, 5, 10, 15),t=3, r=200, p=15, 
                                    tr.method="unif", abun.method="lnorm", 
                                    w.abun=T)

div_sim_res <- as.data.frame(diversity_simulations$results)

diversity_sims_data_frame <- data.frame(div_sim_res$nb.sp,div_sim_res$FDis,
                                        div_sim_res$FDiv,div_sim_res$FEve,
                                        div_sim_res$FRic)
names(diversity_sims_data_frame) <- c("SpeciesNumber","FDis","FDiv","FEve","FRic")
pairs(diversity_sims_data_frame[,-1], upper.panel=NULL)

par(mfrow=c(2,2))
boxplot(FDis~SpeciesNumber,data=diversity_sims_data_frame,ylab="Functional dispersion (FDis)")
boxplot(FDiv~SpeciesNumber,data=diversity_sims_data_frame,ylab="Functional divergence (FDiv)")
boxplot(FEve~SpeciesNumber,data=diversity_sims_data_frame,ylab="Functional evenness (FEve)")
boxplot(FRic~SpeciesNumber,data=diversity_sims_data_frame,ylab="Functional richness (FRic)")
