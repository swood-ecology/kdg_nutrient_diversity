library(haven)                                        # For working with Stata data
library(RStata)                                       # For working with Stata data
library(FD)                                           # For calculating FD metrics
library(stargazer)                                    # To export html code for regression tables
source("~/Box Sync/Work/Statistics/Code/setRows.R")   # for setting rows

# Production data
RICARDIANDATA <- read_dta("~/Box Sync/Work/Writing/Manuscripts/Published/Wood (J Appl Ecol 2017)/RICARDIANDATA.dta")

# Nutrient data
nutr <- read.csv("~/Box Sync/Work/Writing/Manuscripts/Published/Wood (J Appl Ecol 2017)/nutrientData.csv")

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

# Compare metrics to PNA
pna <- read.csv("~/Box Sync/Work/Writing/Manuscripts/Published/Wood (J Appl Ecol 2017)/J Appl Ecol Submission/PNA Analyses/PNA.csv")
diversity_data_frame <- cbind(diversity_data_frame,pna)

# Statistical tests based on ANOVA
# Create regression table
fdis.model <- lm(FDis~SpeciesNumber,data=diversity_data_frame)
fdiv.model <- lm(FDiv~SpeciesNumber,data=diversity_data_frame)
feve.model <- lm(FEve~SpeciesNumber,data=diversity_data_frame)
fric.model <- lm(FRic~SpeciesNumber,data=diversity_data_frame)
pna.model <- lm(PNA~SpeciesNumber,data=diversity_data_frame)

stargazer(fdis.model,fdiv.model,feve.model,fric.model,pna.model,type="html",
          style="aer",omit.stat="f",
          dep.var.labels = c("Functional Dispersion","Functional Divergence",
                             "Functional Evenness","Functional Richness",
                             "Potential Nutrient Adequacy"),
          covariate.labels = c("Two crops","Three crops","Four crops"),
          out="~/Box Sync/Work/Writing/Manuscripts/Submitted/Fonio nutritional diversity/J Appl Ecol Submission/RegTable.html")

# Make plots
# Create function to determine correlations to plot on upper panels
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  if(missing(cex.cor)) cex.cor <- 1/(strwidth(txt))
  cex.final = cex.cor * r
  if(cex.final < .5) cex.final <- .6
  text(0.5, 0.5, txt, cex = cex.final)
}

# Plot diversity against species number
par(mfrow=c(3,2))
boxplot(FDis~SpeciesNumber,data=diversity_data_frame,ylab="Functional dispersion (FDis)")
boxplot(FDiv~SpeciesNumber,data=diversity_data_frame,ylab="Functional divergence (FDiv)")
boxplot(FEve~SpeciesNumber,data=diversity_data_frame,ylab="Functional evenness (FEve)")
boxplot(FRic~SpeciesNumber,data=diversity_data_frame,ylab="Functional richness (FRic)")
boxplot(PNA~SpeciesNumber,data=diversity_data_frame,ylab="Potential Nutrient Adequacy (PNA)")
par(mfrow=c(1,1))


# Plot correlations
pairs(diversity_data_frame[,-c(1,6,8:12)], upper.panel = panel.cor)


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


div_sim_res_low <- as.data.frame(diversity_simulations_low_traits$results)

diversity_sims_low_data_frame <- data.frame(div_sim_res_low$nb.sp,div_sim_res$FDis,
                                        div_sim_res_low$FDiv,div_sim_res$FEve,
                                        div_sim_res_low$FRic)
names(diversity_sims_low_data_frame) <- c("SpeciesNumber","FDis","FDiv","FEve","FRic")
pairs(diversity_sims_low_data_frame[,-1], upper.panel=NULL)

par(mfrow=c(2,2))
boxplot(FDis~SpeciesNumber,data=diversity_sims_low_data_frame,ylab="Functional dispersion (FDis)")
boxplot(FDiv~SpeciesNumber,data=diversity_sims_low_data_frame,ylab="Functional divergence (FDiv)")
boxplot(FEve~SpeciesNumber,data=diversity_sims_low_data_frame,ylab="Functional evenness (FEve)")
boxplot(FRic~SpeciesNumber,data=diversity_sims_low_data_frame,ylab="Functional richness (FRic)")
