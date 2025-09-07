rm(list = ls())
library(tidyverse)
library(multcomp)
library(hommel)
library(matrixStats)
library(mirt)
library(latticeExtra)
library(GMCM) # to use cummean: cumulative mean values


################################################################################
# Parameter specifications for performing statistical tests
################################################################################

nitem <- 10 # Number of items
alphaa <- 0.025 # Nominal significance level

################################################################################
# Reading the PSP dataset
################################################################################
getwd()
list.files()


datPSP <- read_csv("H:/2023_08_PSP_DD/PSP_updated.csv")


### filter subset of the dataset which is of interest (10-item FDA suggested subset of the PSPRS scale from the AbbVie dataset)
datAbbVie <- datPSP %>% filter(STUDY_ID==1,TIME>=0,SCALE==1,VISIT != 21,ITEM %in% c(3, 4, 5, 12, 13, 24, 25, 26, 27, 28))
datAbbVie <- datAbbVie %>% mutate(MDV= ifelse(ITEM>28 | is.na(DV), 1, 0),TIME= ifelse(VISIT==1, 0, TIME))

datAbbVie <- datAbbVie[order(datAbbVie$STUDY_ID, datAbbVie$ID, datAbbVie$ITEM, datAbbVie$VISIT),]

datAbbVie <- datAbbVie %>% mutate(ITEM = case_when(ITEM==3 ~ 1,
                                                   ITEM==4 ~ 2,
                                                   ITEM==5 ~ 3,
                                                   ITEM==12 ~ 4,
                                                   ITEM==13 ~ 5,
                                                   ITEM==24 ~ 6,
                                                   ITEM==25 ~ 7,
                                                   ITEM==26 ~ 8,
                                                   ITEM==27 ~ 9,
                                                   ITEM==28 ~ 10,
                                                   TRUE ~ ITEM),
                                  RICH= ifelse(PHENOTYPE==9, 1, 0), 
                                  INTERV= ifelse(STUDY_ID %in% c(4, 5), 0, 1), 
                                  AGEonset= AGE- DD)



################################################################################
################################################################################
# Estimation of the IRT and linear models for the aggregated data of ABBV-8E12 
# trial (This part of the code is the same in all simulation studies)
################################################################################

### changing the dataset, with visits of interest 1, 6, 11, 14 and 18 (which stands for the baseline and weeks 12, 24, 36 and 52, respectively), from the large to the wide format

datAbbVieALL.wide <- datAbbVie %>% filter(VISIT %in% c(1,6,11,14,18)) %>% 
  pivot_wider(id_cols = c(ID,TRT), names_from = c(ITEM,VISIT), values_from = DV, names_sep = ".", names_prefix = "ITEMV")



### select the data at baseline (only complete cases where the data for all items is available) so that residuals and fitted values can be evaluated in further steps
datAbbVie1.wide <- datAbbVieALL.wide %>% dplyr::select(ends_with(".1")) 
datAbbVie1.wideC <- datAbbVie1.wide[complete.cases(datAbbVie1.wide),]  #361  10

colnames(datAbbVie1.wideC) <- paste0("ITEM",1:10)



### select the data at visit 6 (only complete cases where the data for all items is available), similarly
datAbbVie6.wide <- datAbbVieALL.wide %>% dplyr::select(ends_with(".6")) 
datAbbVie6.wideC <- datAbbVie6.wide[complete.cases(datAbbVie6.wide),]  #360  10

colnames(datAbbVie6.wideC) <- paste0("ITEM",1:10)



### select the data at visit 11 (only complete cases where the data for all items is available), similarly
datAbbVie11.wide <- datAbbVieALL.wide %>% dplyr::select(ends_with(".11")) 
datAbbVie11.wideC <- datAbbVie11.wide[complete.cases(datAbbVie11.wide),]  #350  10

colnames(datAbbVie11.wideC) <- paste0("ITEM",1:10)



### select the data at visit 14 (only complete cases where the data for all items is available), similarly
datAbbVie14.wide <- datAbbVieALL.wide %>% dplyr::select(ends_with(".14")) 
datAbbVie14.wideC <- datAbbVie14.wide[complete.cases(datAbbVie14.wide),]  #243  10

colnames(datAbbVie14.wideC) <- paste0("ITEM",1:10)



### select the data at visit 18 (only complete cases where the data for all items is available), similarly
datAbbVie18.wide <- datAbbVieALL.wide %>% dplyr::select(ends_with(".18")) 
datAbbVie18.wideC <- datAbbVie18.wide[complete.cases(datAbbVie18.wide),] #212  10

colnames(datAbbVie18.wideC) <- paste0("ITEM",1:10)



### combining/pooling the selected data from visits of interest to enhance the precision of the IRT model estimates
datAbbViepull.wideCor <- rbind(datAbbVie1.wideC,datAbbVie6.wideC,datAbbVie11.wideC,datAbbVie14.wideC,datAbbVie18.wideC) #1526   10



### Observe the unique values of each item before rescoring
apply(datAbbViepull.wideCor, 2, function(x) unique(x))



### Rescoring the items based on the FDA recommendation. Note that items 3 and 13, recoded as 1 and 5 in implementations, don't contain any changes.

datAbbViepull.wideC <- datAbbViepull.wideCor %>%
  mutate(
    ITEM2 = case_when(
      ITEM2 == 0 ~ 0,
      ITEM2 == 1 ~ 1,
      ITEM2 == 2 ~ 2,
      ITEM2 == 3 ~ 3,
      ITEM2 == 4 ~ 3,
      TRUE ~ ITEM2
    ),
    ITEM3 = case_when(
      ITEM3 == 0 ~ 0,
      ITEM3 == 1 ~ 1,
      ITEM3 == 2 ~ 1,
      ITEM3 == 3 ~ 1,
      ITEM3 == 4 ~ 2,
      TRUE ~ ITEM3
    ),
    ITEM4 = case_when(
      ITEM4 == 0 ~ 0,
      ITEM4 == 1 ~ 1,
      ITEM4 == 2 ~ 1,
      ITEM4 == 3 ~ 2,
      ITEM4 == 4 ~ 2,
      TRUE ~ ITEM4
    ),
    ITEM6 = case_when(
      ITEM6 == 0 ~ 0,
      ITEM6 == 1 ~ 1,
      ITEM6 == 2 ~ 1,
      ITEM6 == 3 ~ 2,
      ITEM6 == 4 ~ 3,
      TRUE ~ ITEM6
    ),
    ITEM7 = case_when(
      ITEM7 == 0 ~ 0,
      ITEM7 == 1 ~ 0,
      ITEM7 == 2 ~ 0,
      ITEM7 == 3 ~ 1,
      ITEM7 == 4 ~ 2,
      TRUE ~ ITEM7
    ),
    ITEM8 = case_when(
      ITEM8 == 0 ~ 0,
      ITEM8 == 1 ~ 0,
      ITEM8 == 2 ~ 1,
      ITEM8 == 3 ~ 2,
      ITEM8 == 4 ~ 2,
      TRUE ~ ITEM8
    ),
    ITEM9 = case_when(
      ITEM9 == 0 ~ 0,
      ITEM9 == 1 ~ 0,
      ITEM9 == 2 ~ 1,
      ITEM9 == 3 ~ 2,
      ITEM9 == 4 ~ 3,
      TRUE ~ ITEM9
    ),
    ITEM10 = case_when(
      ITEM10 == 0 ~ 0,
      ITEM10 == 1 ~ 0,
      ITEM10 == 2 ~ 1,
      ITEM10 == 3 ~ 2,
      ITEM10 == 4 ~ 3,
      TRUE ~ ITEM10
    ) )


################################################################################

### Observe the unique values of each item after rescoring
apply(datAbbViepull.wideC, 2, unique)



### Fit the graded response IRT model to the dataset of the ABBV-8E12 trial, where aggregated data with original scores from all treatment groups and visits are used for estimation

num.categpull <- as.vector(apply(datAbbViepull.wideCor, 2, function(x) length(unique(x)))) # get the number of categories for each item

IRTmodelAbb.pullC <- mirt(data=datAbbViepull.wideCor,model=1, itemtype = "graded",GR = num.categpull)

res.pullC <- coef(IRTmodelAbb.pullC,simplify=T,IRTpars=T) # estimation of the IRT model parameters

###

### Approximation of the latent variable based on a linear model fit, where as above the aggregated data from ABBV-8E12 trial are considered
### The dependent variable, estimated latent variable from the IRT fitted model, is transformed using the standard logistic function to improve model fit

est.latentAbb.pullCpf <- plogis(fscores(IRTmodelAbb.pullC))
datAbbViepull.wideCor.PSI <- datAbbViepull.wideCor %>% mutate(est.latentAbb.pullCpf=est.latentAbb.pullCpf)

Linmodelpull.CPSIpf <- lm(est.latentAbb.pullCpf ~ ITEM1 + ITEM2 + ITEM3 + ITEM4 +ITEM5 + ITEM6 +ITEM7 + 
                            ITEM8 +ITEM9 + ITEM10, data = datAbbViepull.wideCor.PSI)



### Scatterplot of the estimated latent traits from the IRT and the linear model fit (Left panel of Figure 4 of the paper)

plot(est.latentAbb.pullCpf,fitted(Linmodelpull.CPSIpf), xlab = "Latent trait estimate (transformed)", ylab = "Fitted values",
     main = "Fitted vs. Latent trait (transformed), original score",xlim = c(0, 1), ylim = c(0, 1))

################################################################################

### Similar estimation of the graded response IRT model and the linear model fit
### for the rescored aggregated data from ABBV-8E12 trial

num.categpull.RESC <- as.vector(apply(datAbbViepull.wideC, 2, function(x) length(unique(x)))) # get the number of categories for each item

IRTmodelAbb.pullC.RESC <- mirt(data=datAbbViepull.wideC,model=1, itemtype = "graded",GR = num.categpull.RESC)

res.pullC.RESC <- coef(IRTmodelAbb.pullC.RESC,simplify=T,IRTpars=T) # estimation of the IRT model parameters


### Approximation of the latent variable based on a linear model fit
est.latentAbb.pullCpf.RESC <- plogis(fscores(IRTmodelAbb.pullC.RESC))
datAbbViepull.wideC.PSI <- datAbbViepull.wideC %>% mutate(est.latentAbb.pullCpf.RESC=est.latentAbb.pullCpf.RESC)


Linmodelpull.CPSIpf.RESC <- lm(est.latentAbb.pullCpf.RESC ~ ITEM1 + ITEM2 + ITEM3 + ITEM4 +ITEM5 + ITEM6 +ITEM7 + 
                                 ITEM8 +ITEM9 + ITEM10, data = datAbbViepull.wideC.PSI)



### Scatterplot of the estimated latent traits from the IRT and the linear model fit (Right panel of Figure 4 of the paper)

plot(est.latentAbb.pullCpf.RESC,fitted(Linmodelpull.CPSIpf.RESC), xlab = "Latent trait estimate (transformed)",
     ylab = "Fitted values", main = "Fitted vs. Latent trait (transformed), FDA rescore",xlim = c(0, 1), ylim = c(0, 1))



################################################################################

item.names <- c("Dysp.FS", "Use.KF", "Fall", "Dysa", "Dysp.",
                "Neck.Ri", "Ari.FC", "Gait", "Pos.St", "Sit")


### IRT model estimation of parameters based on both original and the FDA item score (Table 1 of Supplementary Material A)
IRT.model.pars <- round(cbind(res.pullC$items,res.pullC.RESC$items),digits = 3)
row.names(IRT.model.pars) <- item.names
IRT.model.pars
################################################################################

### Linear model fit of the latent variable. The plot shows the normalised weights (Figure 3 in the manuscript)


norm.coef <- coef(Linmodelpull.CPSIpf)[-1]/sum(coef(Linmodelpull.CPSIpf)[-1])

df.coef <- data.frame(Item = factor(item.names, levels = item.names), Value = norm.coef)



norm.coef.RESC <- coef(Linmodelpull.CPSIpf.RESC)[-1]/sum(coef(Linmodelpull.CPSIpf.RESC)[-1])

df.coef.RESC <- data.frame(Item = factor(item.names, levels = item.names), Value = norm.coef.RESC)


df.coef$Score <- "Original"
df.coef.RESC$Score <- "Rescore"


df.coef.both <- rbind(df.coef,df.coef.RESC)


item_order <- unique(df.coef.both$Item)


ggplot(df.coef.both, aes(x = Item, y = Value, color = as.factor(Score))) +
  geom_point(size=2) +
  labs(x = "Items", y = "Norm. estimated coefficients") +
  theme_minimal() +
  scale_color_discrete(name = "Score") +
  scale_x_discrete(limits = item_order) 

