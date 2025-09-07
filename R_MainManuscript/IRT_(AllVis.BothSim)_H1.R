rm(list = ls())
library(tidyverse)
library(multcomp)
library(hommel)
library(matrixStats)
library(mirt)
library(latticeExtra)
library(xtable) #to change tables from R to latex
library(GMCM) # to use cummean: cumulative mean values

library(future)
library(future.apply)  # For parallelized 'apply' functions
library(purrr) 
library(furrr)
library(parallel)


################################################################################
# Parameter specifications for simulation studies
################################################################################

n1 <- 70 # Total sample size for the treatment group
n2 <- 70 # Total sample size for the control group
nitem <- 10 # Number of items 
nmean <- 70 # Total sample per treatment
alphaa <- 0.025 # Nominal significance level
n.sim1 <- 10000 # Number of simulation runs for estimation of the type one error rate
n.sim2 <- 10000 # Number of simulation runs for estimation of the power


################################################################################
# Null distribution of the p-value for the omnibus test
################################################################################

csStat.H0.func <- function(m,N.sim) 
{
  stat.H0.p <- matrix(runif(N.sim*m), nrow = m)
  stat.H0 <- 1/stat.H0.p # the data is (reciprocal) transformed
  stat.H0sort <- apply(stat.H0,2,sort,decreasing=TRUE) #sorting
  csStat.H0 <- apply(stat.H0sort,2,cummean) # cumulative mean values are calculated
  csStat.H0
}

# null distribution of H0 for 10 items
set.seed(2404)
csStat.H0 <-csStat.H0.func(m=10,N.sim=100000)

# null distribution of H0 for domain scores, relating to the sum-score test in three domains of the FDA-10 item version of the PSPRS
set.seed(2310)
csStat.H03dom <-csStat.H0.func(m=3,N.sim=100000)


################################################################################
# Omnibus test as a function of p-values and the null distribution,
# (in principle global, this is used for closed test)
################################################################################

omnibus.test  <- function (p, csStat.H00)
{
  m<-length(p)
  
  p <- matrix(p, ncol = 1)
  stat <- 1/p
  
  csStat <- matrixStats::colCumsums(apply(stat, 2, sort, decreasing = TRUE))/1:m #cummeans
  for (h in 1:m) {
    csStat[h, ] <- ecdf(csStat.H00[h, ])(csStat[h, ]) #rank of observed cumsums
    csStat.H00[h, ] <- rank(csStat.H00[h, ])/ncol(csStat.H00) #rank of the observed cumsums unter H0
  }
  1 - ecdf(colMaxs(csStat.H00))(colMaxs(csStat)) # ) #colMaxs(csStat) Teststat... rank of obs. cumsums among all ranks under H0
}


################################################################################
# Reading the PSP dataset
################################################################################
getwd()
list.files()


#datPSP <- read_csv("H:/2023_08_PSP_DD/PSP_updated.csv")
datPSP <- read_csv("2023_08_PSP_DD/PSP_updated.csv")

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




### Observe the unique values of each item after rescoring
apply(datAbbViepull.wideC, 2, unique)



### Fit the graded response IRT model to the dataset of the ABBV-8E12 trial, where aggregated data with original scores from all treatment groups and visits are used for estimation
num.categpull <- as.vector(apply(datAbbViepull.wideCor, 2, function(x) length(unique(x)))) # get the number of categories for each item

IRTmodelAbb.pullC <- mirt(data=datAbbViepull.wideCor,model=1, itemtype = "graded",GR = num.categpull) # fitting the IRT model

res.pullC <- coef(IRTmodelAbb.pullC,simplify=T,IRTpars=T) # estimation of the IRT model parameters

###

### Approximation of the latent variable based on a linear model fit, where as above the aggregated data from ABBV-8E12 trial are considered
### The dependent variable, estimated latent variable from the IRT fitted model, is transformed using the standard logistic function to improve model fit
est.latentAbb.pullCpf <- plogis(fscores(IRTmodelAbb.pullC))
datAbbViepull.wideCor.PSI <- datAbbViepull.wideCor %>% mutate(est.latentAbb.pullCpf=est.latentAbb.pullCpf)

Linmodelpull.CPSIpf <- lm(est.latentAbb.pullCpf ~ ITEM1 + ITEM2 + ITEM3 + ITEM4 +ITEM5 + ITEM6 +ITEM7 + 
                            ITEM8 +ITEM9 + ITEM10, data = datAbbViepull.wideCor.PSI)


coef(Linmodelpull.CPSIpf) 
fitted.pullpf <- predict(Linmodelpull.CPSIpf, newdata = datAbbViepull.wideCor)


######### Similar estimation of the graded response IRT model and the linear model fit for the rescored aggregated data from ABBV-8E12 trial

num.categpull.RESC <- as.vector(apply(datAbbViepull.wideC, 2, function(x) length(unique(x)))) # get the number of categories for each item

IRTmodelAbb.pullC.RESC <- mirt(data=datAbbViepull.wideC,model=1, itemtype = "graded",GR = num.categpull.RESC) # fitting the IRT model

res.pullC.RESC <- coef(IRTmodelAbb.pullC.RESC,simplify=T,IRTpars=T) # estimation of the IRT model parameters

###

est.latentAbb.pullCpf.RESC <- plogis(fscores(IRTmodelAbb.pullC.RESC))
datAbbViepull.wideC.PSI <- datAbbViepull.wideC %>% mutate(est.latentAbb.pullCpf.RESC=est.latentAbb.pullCpf.RESC)


Linmodelpull.CPSIpf.RESC <- lm(est.latentAbb.pullCpf.RESC ~ ITEM1 + ITEM2 + ITEM3 + ITEM4 +ITEM5 + ITEM6 +ITEM7 + 
                           ITEM8 +ITEM9 + ITEM10, data = datAbbViepull.wideC.PSI)


coef(Linmodelpull.CPSIpf)
fitted.pullpf.RESC <- predict(Linmodelpull.CPSIpf.RESC, newdata = datAbbViepull.wideC)


#############################################################
# Derive relevant Patient IDs for the IRT based simulation
#############################################################

EFFECTRatio <- seq(0.45,0.75,by=0.05)
set.nameEFF <- setNames(EFFECTRatio,paste0("d",1:length(EFFECTRatio)))

### Unique values of patient IDs of the ABBV-8E12 trial dataset at visit 18 (week 52)
ID.summary <- datAbbVie %>% group_by(ID) %>% filter(VISIT==18) %>% pull(ID) %>% unique()

### Unique values of patient IDs of the ABBV-8E12 trial dataset at baseline (week 52)
ID.summary2 <- datAbbVie %>% group_by(ID) %>% filter(VISIT==1) %>% pull(ID) %>% unique()

### intersection of the patient IDs from the two visits. IRT based simulations are based on sampling from these patient IDs and then using their corresponding data
ID.overlap <- intersect(ID.summary,ID.summary2)


###########################################################

# The function which first simulates the IRT based data
# and then performs all the statistical tests for the
# data with Original and rescored values

###########################################################

### This is the function which first simulates the IRT based item scores and then performs all the statistical tests mentioned in Table 2 of the paper
### The function receives all the previous functions and model fits in this script as input for parallel calculation
### the required packages are loaded as well for parallel calculation

TEST.IRIfun <- function(ni,nj,dat, effectRatio = NULL,
                        csStat.H0=csStat.H0,csStat.H03dom=csStat.H03dom,
                        IRTmodelAbb.pullC = IRTmodelAbb.pullC,
                        Linmodelpull.CPSIf = Linmodelpull.CPSIf,
                        Linmodelpull.CPSIpf = Linmodelpull.CPSIpf,
                        IRTmodelAbb.pullC.RESC = IRTmodelAbb.pullC.RESC,
                        Linmodelpull.CPSIf.RESC = Linmodelpull.CPSIf.RESC,
                        Linmodelpull.CPSIpf.RESC = Linmodelpull.CPSIpf.RESC){
  
  
  ### sampling ni+nj=140 patient IDs from the "ID.overlap" defined above 
  v<- sample(ID.overlap,size = ni+nj,replace = FALSE)
  dat2 <- dat %>% dplyr::filter(ID %in% v)
  
  
  ### Setting new treatment labels for the test and control groups in the sampled data, as the original labels contain all three groups 1,2,3
  TRT.lable <- sample(c(rep(1,ni),rep(2,nj)),size = ni+nj,replace = FALSE)
  val.names <- setNames(TRT.lable,unique(dat2$ID))
  dat2 <- dat2 %>%  mutate(TRT.NEW= val.names[match(ID,names(val.names))])
  
  dat2 <- dat2 %>% dplyr::select(-c(NID, FD)) 
  dat2 <- dat2 %>% mutate(TIME=TIME/365) #time variable in year
  ###########################################
  
  PSI_MODEL=0
  UNDEF=0
  OC1=1   
  
  # set seed here again for reproducibilty
  # set.seed()
  dat2 <- dat2 %>% dplyr::group_by(ID) %>% dplyr::mutate(ETA1=rnorm(1, 0,sqrt(0.666888)), ETA2= rnorm(1,0, sqrt(0.365698)))
  
  ###########################################
  ### A two-parameter IRT model is fitted to the PSPRS FDA subitems and then the Item characteristic parameters (i.e., discrimination and difficulty)
  ### for each of the 10 items are estimated. Those are the fixed effects parameters THETA1 to THETA50. The estimation is done using NONMEM software.
  
  THETA1 <- 1.0551 # ; I1DIS 1
  THETA2 <- -0.51518 # ; I1DIF1 2
  THETA3 <- 2.759 # ; I1DIF2 3
  THETA4 <- 1.4513 # ; I1DIF3 4
  THETA5 <- 1.8269 # ; I1DIF4 5
  THETA6 <- 1.7171 # ; I2DIS 6
  THETA7 <- -2.031 # ; I2DIF1 7
  THETA8 <- 1.6387 # ; I2DIF2 8
  THETA9 <- 1.0799 # ; I2DIF3 9
  THETA10 <- 1.7038 # ; I2DIF4 10
  THETA11 <- 1.1596 # ; I3DIS 11
  THETA12 <- -2.7056 # ; I3DIF1 12
  THETA13 <- 1.4857 # ; I3DIF2 13
  THETA14 <- 1.3179 # ; I3DIF3 14
  THETA15 <- 1.2065 # ; I3DIF4 15
  THETA16 <- 1.1677 # ; I4DIS 16
  THETA17 <- -2.6012 # ; I4DIF1 17
  THETA18 <- 2.2392 # ; I4DIF2 18
  THETA19 <- 1.9975 # ; I4DIF3 19
  THETA20 <- 1.5462 # ; I4DIF4 20
  THETA21 <- 0.98833 # ; I5DIS 21
  THETA22 <- -1.1392 # ; I5DIF1 22
  THETA23 <- 1.4591 # ; I5DIF2 23
  THETA24 <- 1.593 # ; I5DIF3 24
  THETA25 <- 3.4615 # ; I5DIF4 25
  THETA26 <- 1.061 # ; I6DIS 26
  THETA27 <- -2.1605 # ; I6DIF1 27
  THETA28 <- 1.5669 # ; I6DIF2 28
  THETA29 <- 1.7424 # ; I6DIF3 29
  THETA30 <- 2.2099 # ; I6DIF4 30
  THETA31 <- 3.4247 # ; I7DIS 31
  THETA32 <- -1.4049 # ; I7DIF1 32
  THETA33 <- 0.91454 # ; I7DIF2 33
  THETA34 <- 0.32321 # ; I7DIF3 34
  THETA35 <- 0.9472 # ; I7DIF4 35
  THETA36 <- 3.7782 # ; I8DIS 36
  THETA37 <- -2.0797 # ; I8DIF1 37
  THETA38 <- 1.4038 # ; I8DIF2 38
  THETA39 <- 0.83733 # ; I8DIF3 39
  THETA40 <- 1.5922 # ; I8DIF4 40
  THETA41 <- 2.768 # ; I9DIS 41
  THETA42 <- -1.743 # ; I9DIF1 42
  THETA43 <- 0.87472 # ; I9DIF2 43
  THETA44 <- 0.77642 # ; I9DIF3 44
  THETA45 <- 0.96647 # ; I9DIF4 45
  THETA46 <- 3.5292 # ; I10DIS 46
  THETA47 <- -1.4802 # ; I10DIF1 47
  THETA48 <- 1.0276 # ; I10DIF2 48
  THETA49 <- 0.83825 # ; I10DIF3 49
  THETA50 <- 1.1534 # ; I10DIF4 50
  
  
  
  
  ### exploring the covariate model:  a step wise algorithm has been used to find the combination of covariates that 
  ### result in the most significant improvement in model fit using NONMEM software. It resulted into the covariates THETA53-THETA57
  
  ### A simple longitudinal model has been explored to describe the progression of the latent variable using a single baseline and slope parameters
  ### (THETA51 and THETA52) with interindividual variability (ETA1 and ETA2) 
  
  THETA51 <- -0.311627 #(fixed effect) estimation of the baseline
  THETA52 <- 0.884462 #(fixed effect) estimation of the slope
  THETA53 <- -0.0304938 # BASEAGE1
  THETA54 <- 0.877237 # BASEINTERV1
  THETA55 <- 1.16983 # BASERICH1
  THETA56 <- -0.358592 # BASESEX1
  THETA57 <- -0.264083 # SLOPEINTERV1
  
  
  EPS1 <- 1
  
  ###########################################
  
  dat2 <- dat2 %>% dplyr::mutate(BASESEX=ifelse(SEX==2, 1, 1 + THETA56))
  
  
  dat2 <- dat2 %>% dplyr::mutate(SLOPEINTERV=ifelse(INTERV==1, 1, 1 + THETA57))
  
  
  dat2 <- dat2 %>% dplyr::mutate(SLOPECOV=ifelse(is.na(SLOPEINTERV),1, SLOPEINTERV))
  
  
  dat2 <- dat2 %>% dplyr::mutate(BASERICH=ifelse(RICH==1, 1, 1 + THETA55))
  
  
  dat2 <- dat2 %>% dplyr::mutate(BASEINTERV=ifelse(INTERV==1, 1, 1 + THETA54))
  
  
  dat2 <- dat2 %>% dplyr::mutate(BASEAGE=( 1 + THETA53*(AGE - 69)))
  
  
  
  ### Assignment of item parameters
  
  dat2 <- dat2 %>% dplyr::ungroup()
  dat2 <- dat2 %>% dplyr::mutate(DIS = case_when(ITEM==1 ~ THETA1,
                                                 ITEM==2 ~ THETA6,
                                                 ITEM==3 ~ THETA11,
                                                 ITEM==4 ~ THETA16,
                                                 ITEM==5 ~ THETA21,
                                                 ITEM==6 ~ THETA26,
                                                 ITEM==7 ~ THETA31,
                                                 ITEM==8 ~ THETA36,
                                                 ITEM==9 ~ THETA41,
                                                 ITEM==10 ~ THETA46,
                                                 TRUE ~ 0))
  
  
  
  dat2 <- dat2 %>% dplyr::mutate(DIF1 = case_when(ITEM==1 ~ THETA2,
                                                  ITEM==2 ~ THETA7,
                                                  ITEM==3 ~ THETA12,
                                                  ITEM==4 ~ THETA17,
                                                  ITEM==5 ~ THETA22,
                                                  ITEM==6 ~ THETA27,
                                                  ITEM==7 ~ THETA32,
                                                  ITEM==8 ~ THETA37,
                                                  ITEM==9 ~ THETA42,
                                                  ITEM==10 ~ THETA47,
                                                  TRUE ~ 0))
  
  
  dat2 <- dat2 %>% dplyr::mutate(DIF2 = case_when(ITEM==1 ~ THETA3,
                                                  ITEM==2 ~ THETA8,
                                                  ITEM==3 ~ THETA13,
                                                  ITEM==4 ~ THETA18,
                                                  ITEM==5 ~ THETA23,
                                                  ITEM==6 ~ THETA28,
                                                  ITEM==7 ~ THETA33,
                                                  ITEM==8 ~ THETA38,
                                                  ITEM==9 ~ THETA43,
                                                  ITEM==10 ~ THETA48,
                                                  TRUE ~ 0))
  
  
  dat2 <- dat2 %>% dplyr::mutate(DIF3 = case_when(ITEM==1 ~ THETA4,
                                                  ITEM==2 ~ THETA9,
                                                  ITEM==3 ~ THETA14,
                                                  ITEM==4 ~ THETA19,
                                                  ITEM==5 ~ THETA24,
                                                  ITEM==6 ~ THETA29,
                                                  ITEM==7 ~ THETA34,
                                                  ITEM==8 ~ THETA39,
                                                  ITEM==9 ~ THETA44,
                                                  ITEM==10 ~ THETA49,
                                                  TRUE ~ 0))
  
  
  dat2 <- dat2 %>% dplyr::mutate(DIF4 = case_when(ITEM==1 ~ THETA5,
                                                  ITEM==2 ~ THETA10,
                                                  ITEM==3 ~ THETA15,
                                                  ITEM==4 ~ THETA20,
                                                  ITEM==5 ~ THETA25,
                                                  ITEM==6 ~ THETA30,
                                                  ITEM==7 ~ THETA35,
                                                  ITEM==8 ~ THETA40,
                                                  ITEM==9 ~ THETA45,
                                                  ITEM==10 ~ THETA50,
                                                  TRUE ~ 0))
  
  
  
  ###########################################
  
  ### Constructing the baseline and slope parameters in the longitudinal model
  
  dat2 <- dat2 %>% dplyr::mutate(BASECOV=BASEAGE*BASEINTERV*BASERICH*BASESEX)
  dat2 <- dat2 %>% dplyr::mutate(BASECOV=is.na(BASECOV), 1, BASECOV)
  
  dat2 <- dat2 %>% dplyr::mutate(TVBASE= THETA51)
  
  dat2 <- dat2 %>% dplyr::mutate(TVBASE = BASECOV*TVBASE)
  
  dat2 <- dat2 %>% dplyr::mutate(BASE= TVBASE+ETA1)
  
  dat2 <- dat2 %>% dplyr::mutate(TVSLOPE=THETA52)
  
  dat2 <- dat2 %>% dplyr::mutate(TVSLOPE = SLOPECOV*TVSLOPE)
  
  dat2 <- dat2 %>% dplyr::mutate(SLP=TVSLOPE+ETA2)
  
  
  ### The progression model of the latent variable
  
  if(is.null(effectRatio)){
    
    dat2 <- dat2 %>% dplyr::mutate(PSI = BASE + SLP*TIME)
    
  }else{
    
    dat2 <- dat2 %>% dplyr::mutate(PSI = ifelse(TRT.NEW==2 & VISIT!=1 , BASE + effectRatio*SLP*TIME , BASE + SLP*TIME) )
    
  }
  
  ### Ordered categorical data model with 5 levels: [0,4]
  
  dat2 <- dat2 %>% dplyr::mutate(DIFG1=DIF1)
  dat2 <- dat2 %>% dplyr::mutate(DIFG2=DIFG1+DIF2)
  dat2 <- dat2 %>% dplyr::mutate(DIFG3=DIFG2+DIF3)
  dat2 <- dat2 %>% dplyr::mutate(DIFG4=DIFG3+DIF4)
  
  dat2 <- dat2 %>% dplyr::mutate(PGE1=exp(DIS*(PSI-DIFG1))/(1+exp(DIS*(PSI-DIFG1))))
  dat2 <- dat2 %>% dplyr::mutate(PGE2=exp(DIS*(PSI-DIFG2))/(1+exp(DIS*(PSI-DIFG2))))
  dat2 <- dat2 %>% dplyr::mutate(PGE3=exp(DIS*(PSI-DIFG3))/(1+exp(DIS*(PSI-DIFG3))))
  dat2 <- dat2 %>% dplyr::mutate(PGE4=exp(DIS*(PSI-DIFG4))/(1+exp(DIS*(PSI-DIFG4))))
  
  

  ### Simulation code of data generation using IRT
  
  dat2 <- dat2 %>% dplyr::ungroup() %>% dplyr::mutate(p_sim=runif(nrow(dat2),0,1))
  dat2 <- dat2 %>% dplyr::mutate(DV=0)
  dat2 <- dat2 %>% dplyr::mutate(DV = case_when(p_sim<PGE4 ~ 4,
                                                p_sim<PGE3 ~ 3,
                                                p_sim<PGE2 ~ 2,
                                                p_sim<PGE1 ~ 1,
                                                TRUE ~ DV))
  
  
  
  ### Use the simulated data from two time points (Baseline and week 52) to fit the ANCOVA models
  ### ANCOVA models are fitted to simulated data with the original scores and the FDA rescores
  
  dat2.118 <- dat2 %>% filter(VISIT %in% c(1,18)) 
  
  
  ### Rescoring the items based on the FDA recommendation. Note that items 3 and 13, recoded as 1 and 5 in implementations, don't contain any changes.
  dat2.118rec <- dat2.118 %>%
    mutate(
      DV = case_when( 
        ITEM==1& DV == 0 ~ 0,
        ITEM==1& DV == 1 ~ 1,
        ITEM==1& DV == 2 ~ 2,
        ITEM==1& DV == 3 ~ 3,
        ITEM==1& DV == 4 ~ 4,
        TRUE ~ DV
      ),
      DV = case_when(
        ITEM==2& DV == 0 ~ 0,
        ITEM==2& DV == 1 ~ 1,
        ITEM==2& DV == 2 ~ 2,
        ITEM==2& DV == 3 ~ 3,
        ITEM==2& DV == 4 ~ 3,
        TRUE ~ DV
      ),
      DV = case_when(
        ITEM==3& DV == 0 ~ 0,
        ITEM==3& DV == 1 ~ 1,
        ITEM==3& DV == 2 ~ 1,
        ITEM==3& DV == 3 ~ 1,
        ITEM==3& DV == 4 ~ 2,
        TRUE ~ DV
      ),
      DV = case_when(
        ITEM==4& DV == 0 ~ 0,
        ITEM==4& DV == 1 ~ 1,
        ITEM==4& DV == 2 ~ 1,
        ITEM==4& DV == 3 ~ 2,
        ITEM==4& DV == 4 ~ 2,
        TRUE ~ DV
      ),
      DV = case_when(
        ITEM==5& DV == 0 ~ 0,
        ITEM==5& DV == 1 ~ 1,
        ITEM==5& DV == 2 ~ 2,
        ITEM==5& DV == 3 ~ 3,
        ITEM==5& DV == 4 ~ 4,
        TRUE ~ DV
      ),
      DV = case_when(
        ITEM==6& DV == 0 ~ 0,
        ITEM==6& DV == 1 ~ 1,
        ITEM==6& DV == 2 ~ 1,
        ITEM==6& DV == 3 ~ 2,
        ITEM==6& DV == 4 ~ 3,
        TRUE ~ DV
      ),
      DV = case_when(
        ITEM==7& DV == 0 ~ 0,
        ITEM==7& DV == 1 ~ 0,
        ITEM==7& DV == 2 ~ 0,
        ITEM==7& DV == 3 ~ 1,
        ITEM==7& DV == 4 ~ 2,
        TRUE ~ DV
      ),
      DV = case_when( 
        ITEM==8& DV == 0 ~ 0,
        ITEM==8& DV == 1 ~ 0,
        ITEM==8& DV == 2 ~ 1,
        ITEM==8& DV == 3 ~ 2,
        ITEM==8& DV == 4 ~ 2,
        TRUE ~ DV
      ),
      DV = case_when( 
        ITEM==9& DV == 0 ~ 0,
        ITEM==9& DV == 1 ~ 0,
        ITEM==9& DV == 2 ~ 1,
        ITEM==9& DV == 3 ~ 2,
        ITEM==9& DV == 4 ~ 3,
        TRUE ~ DV
      ),
      DV = case_when( 
        ITEM==10& DV == 0 ~ 0,
        ITEM==10& DV == 1 ~ 0,
        ITEM==10& DV == 2 ~ 1,
        ITEM==10& DV == 3 ~ 2,
        ITEM==10& DV == 4 ~ 3,
        TRUE ~ DV
      ) )
  
  
  
  ### Change the format of simulated data from long to wide to fit ANCOVA models
  
  ### Original data preparation
  dat2.piw <- dat2 %>% filter(VISIT %in% c(1,18)) %>% 
    pivot_wider(id_cols = c(ID,TRT.NEW,SCALE), names_from = c(ITEM,VISIT), values_from = DV, names_prefix = "ITEMV",names_sep = ".") %>%
    relocate(ends_with(".1"), .after = last_col())
  
  dat2.piw$TRT.NEW <- factor(dat2.piw$TRT.NEW,labels = c("1.control","2.test"))
  
  
  TRT.LABLED <- dat2.piw$TRT.NEW
  
  
  
  ### rescored data preparation
  
  dat2.piw.RESC <- dat2.118rec %>% filter(VISIT %in% c(1,18)) %>% 
    pivot_wider(id_cols = c(ID,TRT.NEW,SCALE), names_from = c(ITEM,VISIT), values_from = DV, names_prefix = "ITEMV",names_sep = ".") %>%
    relocate(ends_with(".1"), .after = last_col())
  
  
  
  dat2.piw.RESC <- dat2.piw.RESC %>%
    mutate( TRT.LABLED=TRT.LABLED )
  
  
  ###########################################
  ###########################################
  
  ### Perform all the statistical tests for the simulated data with the original scores
  
  ### Estimation of the latent variable based on the IRT model for the simulated data (with the original scores) at baseline and week 52
  Fscores.est1 <- dat2.piw %>% dplyr::select(ends_with(".1"))  %>%
    apply(., 1, function(row) fscores(IRTmodelAbb.pullC,response.pattern = row)[1])
  
  Fscores.est18 <- dat2.piw %>% dplyr::select(ends_with(".18")) %>%
    apply(., 1, function(row) fscores(IRTmodelAbb.pullC,response.pattern = row)[1])
  
  
  
  ### Fitting the ANCOVA model with the estimated latent variable at week 52 as the dependent variable and
  ### the treatment and the estimated latent trait at baseline as the independent variables
  dat2.piwPSI <- dat2.piw %>%
    mutate(Fscores.est18=Fscores.est18,Fscores.est1=Fscores.est1)
  
  mPSIf <- lm(formula=Fscores.est18 ~ TRT.NEW + Fscores.est1, data = dat2.piwPSI )
  
  
  model.summary <- summary(mPSIf)
  
  # Extract the degrees of freedom of residuals
  df.residual <- model.summary$df[2]
  
  
  ### Approximate IRT-based test
  
  dat2.piw1 <- dat2.piw %>% dplyr::select(ends_with(".1"))
  dat2.piw18 <- dat2.piw %>% dplyr::select(ends_with(".18"))
  
  colnames(dat2.piw1) <- colnames(dat2.piw18) <- paste0("ITEM",1:nitem)
  
  
  ### Approximation of the latent variable based on the linear model fit
  linmPSIpf.est1 <- predict(Linmodelpull.CPSIpf, newdata = dat2.piw1) # from the linear model with input plogis(fscores)
  linmPSIpf.est18 <- predict(Linmodelpull.CPSIpf, newdata = dat2.piw18)
  
  ### Using the standard logit function for back transformation of the predicted latent variables
  linmPSIbpf.est1 <- qlogis(linmPSIpf.est1)
  linmPSIbpf.est18 <- qlogis(linmPSIpf.est18)
  
  dat2.piw.linmPSI <- dat2.piw %>%
    mutate(linmPSIbpf.est1=linmPSIbpf.est1,linmPSIbpf.est18=linmPSIbpf.est18)
  
  
  m.linmPSIbpf <- lm(formula=linmPSIbpf.est18 ~ TRT.NEW + linmPSIbpf.est1, data = dat2.piw.linmPSI ) #back transformed
  
  ####
  
  dat2.piwsum <- dat2.piw %>% mutate(
    
    ### Sum score across items, at baseline and week 52
    sumStest = rowSums(dplyr::select(., ends_with("18"))),
    sumSbase = rowSums(dplyr::select(., ends_with("1"))),
    
    
    ### Sum score across items for HISTORY domain, at baseline and week 52
    sumStest.dom1 = rowSums(dplyr::select(., c(ITEMV1.18,ITEMV2.18,ITEMV3.18))),
    sumSbase.dom1 = rowSums(dplyr::select(., c(ITEMV1.1,ITEMV2.1,ITEMV3.1))),
    
    
    ### Sum score across items for BULBAR EXAM domain, at baseline and week 52
    sumStest.dom2 = rowSums(dplyr::select(., c(ITEMV4.18,ITEMV5.18))),
    sumSbase.dom2 = rowSums(dplyr::select(., c(ITEMV4.1,ITEMV5.1))),
    
    
    ### Sum score across items for GAIT/MIDLINE domain, at baseline and week 52
    sumStest.dom3 = rowSums(dplyr::select(., c(ITEMV6.18,ITEMV7.18,ITEMV8.18,ITEMV9.18,ITEMV10.18 ))),
    sumSbase.dom3 = rowSums(dplyr::select(., c(ITEMV6.1,ITEMV7.1,ITEMV8.1,ITEMV9.1,ITEMV10.1 )))
    
  )
  
  
  # ANCOVA model for the sum scores
  mITEMall <- lm(formula=sumStest ~ TRT.NEW + sumSbase, data = dat2.piwsum )
  
  
  # ANCOVA models for the sum scores across items in each of the three domains
  mITEM.dom1 <- lm(formula=sumStest.dom1 ~ TRT.NEW + sumSbase.dom1, data = dat2.piwsum )
  
  mITEM.dom2 <- lm(formula=sumStest.dom2 ~ TRT.NEW + sumSbase.dom2, data = dat2.piwsum )
  
  mITEM.dom3 <- lm(formula=sumStest.dom3 ~ TRT.NEW + sumSbase.dom3, data = dat2.piwsum )
  
  
  ### P-value of the sum score test for the items in the HISTORY domain 
  pval.SumScore.dom1 <- pt(summary(mITEM.dom1)$coef["TRT.NEW2.test","t value"],df=ni+nj-3,lower.tail = T)
  
  
  ### P-value of the sum score test for the items in the BULBAR EXAM domain 
  pval.SumScore.dom2 <- pt(summary(mITEM.dom2)$coef["TRT.NEW2.test","t value"],df=ni+nj-3,lower.tail = T)
  
  
  ### P-value of the sum score test for the items in the GAIT/MIDLINE domain 
  pval.SumScore.dom3 <- pt(summary(mITEM.dom3)$coef["TRT.NEW2.test","t value"],df=ni+nj-3,lower.tail = T)
  
  
  ### ANCOVA models for individual items (marginal ANCOVA models)
  m3 <-  lm(formula=ITEMV1.18 ~ TRT.NEW + ITEMV1.1, data = dat2.piw )
  m4 <-  lm(formula=ITEMV2.18 ~ TRT.NEW + ITEMV2.1, data = dat2.piw )
  m5 <-  lm(formula=ITEMV3.18 ~ TRT.NEW + ITEMV3.1, data = dat2.piw )
  m12 <- lm(formula=ITEMV4.18 ~TRT.NEW + ITEMV4.1, data = dat2.piw )
  m13 <-  lm(formula=ITEMV5.18 ~ TRT.NEW + ITEMV5.1, data = dat2.piw )
  m24 <-  lm(formula=ITEMV6.18 ~ TRT.NEW + ITEMV6.1, data = dat2.piw )
  m25 <-  lm(formula=ITEMV7.18 ~ TRT.NEW + ITEMV7.1, data = dat2.piw )
  m26 <- lm(formula=ITEMV8.18 ~ TRT.NEW + ITEMV8.1, data = dat2.piw )
  m27 <-  lm(formula=ITEMV9.18 ~ TRT.NEW + ITEMV9.1, data = dat2.piw )
  m28 <-  lm(formula=ITEMV10.18 ~ TRT.NEW + ITEMV10.1, data = dat2.piw )
  
  
  ### Multiple marginal modeling
  MMMall <- mmm(m3,m4,m5,m12,m13,m24,m25,m26,m27,m28)
  
  
  ### A function to define the contrasts for the treatment effect estimation in each marginal ANCOVA model
  contfun <- function(xx) c(rep(0,3*xx-2),1,rep(0,30-(3*xx-2+1)))
  
  mcontf <- matrix( c(contfun(1),contfun(2),contfun(3),contfun(4),contfun(5),
                      contfun(6),contfun(7),contfun(8),contfun(9),contfun(10)), nrow = 10, ncol = 30 ,byrow = TRUE )
  
  
  glhij <- glht(MMMall,linfct = mcontf,alternative="less")
  
  
  ### Estimation of the treatment effect and the covariance function of items from multiple marginal models
  t.eachij <- coef(glhij)/sqrt(diag(vcov(glhij)))
  gencovij.pool <- vcov(glhij)
  gencorij.pool <- cov2cor(gencovij.pool)
  
  
  
  ### P-value for each individual item
  pval.eachij <- numeric(nitem)
  for(k in 1:nitem){
    pval.eachij[k] <- pt(t.eachij[k],df=ni+nj-3,lower.tail = T)
  }
  
  
  ### The OLS test statistic
  t.OLSij <- sum(t.eachij)/sqrt(sum(gencorij.pool))
  
  
  ### Modified approximation of the degrees of freedom for the OLS and GLS tests
  dfOGij <- 0.5*(ni+nj-3)*(1+1/nitem^2)
  
  ### P-value of the OLS test
  pval.OLSij <- pt(t.OLSij,df=dfOGij,lower.tail = T)
  
  ### The weight for each individual item in the OLS test
  WOLSij <- 1/sqrt(sum(gencorij.pool))
  
  
  ### The weight for each individual item in the GLS test
  WGLSij <- colSums(solve(gencorij.pool))/sqrt(sum(solve(gencorij.pool)))
  names(WGLSij) <- c(3,4,5,12,13,24,25,26,27,28)
  
  
  # NegRij records proportions of simulations where at least one negative weight in the GLS test occur (not reported in the paper)
  # negItEM records items with negative weights (not reported in the paper)
  if(length(WGLSij[WGLSij<0])>=1){
    negItEM <- as.integer(names(WGLSij)[WGLSij<0])
    NegRij <- 1
  }else{
    negItEM <- 100
    NegRij <- 0
  }
  
  
  ### The GLS test statistic
  t.GLSij <- sum(solve(gencorij.pool)%*%t.eachij)/sqrt(sum(solve(gencorij.pool)))
  
  
  ### P-value of the GLS test
  pval.GLSij <- pt(t.GLSij,df=dfOGij,lower.tail = T)
  
  
  ### GLS test after dropping item 26 ( Gait), denoted as GLS-26
  gencovij.pool26 <- gencovij.pool[-c(8),-c(8)]
  gencorij.pool26 <- cov2cor(gencovij.pool26)
  
  WGLSij26 <- colSums(solve(gencorij.pool26))/sqrt(sum(solve(gencorij.pool26)))
  names(WGLSij26) <- c(3,4,5,12,13,24,25,27,28)
  
  
  # NegRij26 records proportions of simulations where at least one negative weight in the GLS-26 test occur (not reported in the paper)
  # negItEM26 records items with negative weights in the GLS-26 test (not reported in the paper)
  
  if(length(WGLSij26[WGLSij26<0])>=1){
    negItEM26 <- as.integer(names(WGLSij26)[WGLSij26<0]) 
    NegRij26 <- 1
  }else{
    negItEM26 <- 100
    NegRij26 <- 0
  }
  
  
  t.GLSij26 <- sum(solve(gencorij.pool26)%*%t.eachij[-8])/sqrt(sum(solve(gencorij.pool26)))
  dfOGij26 <- 0.5*(ni+nj-3)*(1+1/(nitem-1)^2)
  pval.GLSij26 <- pt(t.GLSij26,df=dfOGij26,lower.tail = T)
  
  
  ### Rejection rate in the GLS-26 test
  res.GLSij26 <- ifelse(pval.GLSij26<alphaa,1,0)
  
  ####
  
  ### Rejection rates in the OLS and GLS tests
  res.OLSij <- ifelse(pval.OLSij<alphaa,1,0)
  res.GLSij <- ifelse(pval.GLSij<alphaa,1,0)
  
  
  ### Rejection rate in the Bonferroni correction test
  pval.bonfij <- nitem*min(pval.eachij)
  if(pval.bonfij>=1)  pval.bonfij <- 1
  res.bonfij <- ifelse(pval.bonfij<alphaa,1,0)
  
  
  ### P-value and the rejection rate in the Simes test
  pval.Simesij <- min(hommel::p.adjust(hommel(pval.eachij)))
  res.Simesij <- ifelse(pval.Simesij<alphaa,1,0)
  
  
  ### P-value and rejection rate based on minimum T-values (or the corresponding Z-values)
  z.min <- min(qnorm(pt(t.eachij,df=ni+nj-3)))
  pval.tmin <- 1-pmvnorm(lower=rep(z.min,nitem),upper=rep(Inf,nitem),corr = as.matrix(gencorij.pool) )[1]
  res.tminij <- ifelse(pval.tmin<alphaa,1,0)
  
  
  ### P-value and rejection rate in the  Sum Score test
  pval.SumScoreij <- pt(summary(mITEMall)$coef["TRT.NEW2.test","t value"],df=ni+nj-3,lower.tail = T)
  res.SumScoreij <- ifelse(pval.SumScoreij<alphaa,1,0)
  
  
  ### P-value and rejection rate in the IRT-based test (abbreviated as IRT.PSIF)
  pval.PSIfij <- pt(summary(mPSIf)$coef["TRT.NEW2.test","t value"],df=ni+nj-3,lower.tail = T)
  res.PSIfij <- ifelse(pval.PSIfij<alphaa,1,0)
  
  
  ### P-value and rejection rate in the Linear model-based test (abbreviated as LM.PSIBPF)
  pval.linmPSIbpfij <- pt(summary(m.linmPSIbpf)$coef["TRT.NEW2.test","t value"],df=ni+nj-3,lower.tail = T)
  res.linmPSIbpfij <- ifelse(pval.linmPSIbpfij<alphaa,1,0)
  
  
  
  ### P-value and rejection rate in the omnibus test
  pval.Omnibus <- omnibus.test(pval.eachij,csStat.H0)
  res.Omnibusij <- ifelse(pval.Omnibus<alphaa,1,0)
  
  
  
  ### P-value and rejection rate in the combination of the Omnibus and sum score test (abbreviated as Omnibus-dom)
  pval.Omnibus3Domij <- omnibus.test(c(pval.SumScore.dom1,pval.SumScore.dom2,pval.SumScore.dom3),csStat.H03dom)
  res.Omnibus3Domij <- ifelse(pval.Omnibus3Domij<alphaa,1,0)
  
  
  ###########################################
  ###########################################
  
  ### Perform all the statistical tests for the rescored simulated data
  ### All the tests and calculations are the same, only the test names are modified to account for rescoring
  
  Fscores.est1.RESC <- dat2.piw.RESC %>% dplyr::select(ends_with(".1"))  %>%
    apply(., 1, function(row) fscores(IRTmodelAbb.pullC.RESC,response.pattern = row)[1])
  
  Fscores.est18.RESC <- dat2.piw.RESC %>% dplyr::select(ends_with(".18")) %>%
    apply(., 1, function(row) fscores(IRTmodelAbb.pullC.RESC,response.pattern = row)[1])
  
  
  dat2.piwPSI.RESC <- dat2.piw.RESC %>%
    mutate(Fscores.est18.RESC=Fscores.est18.RESC,Fscores.est1.RESC=Fscores.est1.RESC)
  
  mPSIf.RESC <- lm(formula=Fscores.est18.RESC ~ TRT.LABLED + Fscores.est1.RESC, data = dat2.piwPSI.RESC )
  
  
  model.summary.RESC <- summary(mPSIf.RESC)
  
  
  df.residual.RESC <- model.summary.RESC$df[2]
  
  ####
  
  dat2.piw1.RESC <- dat2.piw.RESC %>% dplyr::select(ends_with(".1"))
  dat2.piw18.RESC <- dat2.piw.RESC %>% dplyr::select(ends_with(".18"))
  
  colnames(dat2.piw1.RESC) <- colnames(dat2.piw18.RESC) <- paste0("ITEM",1:nitem)
  
  
  linmPSIpf.est1.RESC <- predict(Linmodelpull.CPSIpf.RESC, newdata = dat2.piw1.RESC)
  linmPSIpf.est18.RESC <- predict(Linmodelpull.CPSIpf.RESC, newdata = dat2.piw18.RESC)
  
  
  linmPSIbpf.est1.RESC <- qlogis(linmPSIpf.est1.RESC)
  linmPSIbpf.est18.RESC <- qlogis(linmPSIpf.est18.RESC)
  
  dat2.piw.linmPSI.RESC <- dat2.piw.RESC %>%
    mutate(linmPSIbpf.est1.RESC=linmPSIbpf.est1.RESC,linmPSIbpf.est18.RESC=linmPSIbpf.est18.RESC,)
  
  
  m.linmPSIbpf.RESC <- lm(formula=linmPSIbpf.est18.RESC ~ TRT.LABLED + linmPSIbpf.est1.RESC, data = dat2.piw.linmPSI.RESC )
  
  #####
  
  dat2.piwsum.RESC <- dat2.piw.RESC %>% mutate(
    sumStest.RESC = rowSums(dplyr::select(., ends_with("18"))),
    sumSbase.RESC = rowSums(dplyr::select(., ends_with("1"))),
    
    sumStest.dom1.RESC = rowSums(dplyr::select(., c(ITEMV1.18,ITEMV2.18,ITEMV3.18))),
    sumSbase.dom1.RESC = rowSums(dplyr::select(., c(ITEMV1.1,ITEMV2.1,ITEMV3.1))),
    
    sumStest.dom2.RESC = rowSums(dplyr::select(., c(ITEMV4.18,ITEMV5.18))),
    sumSbase.dom2.RESC = rowSums(dplyr::select(., c(ITEMV4.1,ITEMV5.1))),
    
    sumStest.dom3.RESC = rowSums(dplyr::select(., c(ITEMV6.18,ITEMV7.18,ITEMV8.18,ITEMV9.18,ITEMV10.18 ))),
    sumSbase.dom3.RESC = rowSums(dplyr::select(., c(ITEMV6.1,ITEMV7.1,ITEMV8.1,ITEMV9.1,ITEMV10.1 )))
    
  )
  

  mITEMall.RESC <- lm(formula=sumStest.RESC ~ TRT.LABLED + sumSbase.RESC, data = dat2.piwsum.RESC )
  
  
  mITEM.dom1.RESC <- lm(formula=sumStest.dom1.RESC ~ TRT.LABLED + sumSbase.dom1.RESC, data = dat2.piwsum.RESC )
  
  
  mITEM.dom2.RESC <- lm(formula=sumStest.dom2.RESC ~ TRT.LABLED + sumSbase.dom2.RESC, data = dat2.piwsum.RESC )
  
  
  mITEM.dom3.RESC <- lm(formula=sumStest.dom3.RESC ~ TRT.LABLED + sumSbase.dom3.RESC, data = dat2.piwsum.RESC )

  
  pval.SumScore.dom1.RESC <- pt(summary(mITEM.dom1.RESC)$coef["TRT.LABLED2.test","t value"],df=ni+nj-3,lower.tail = T)
  
  
  pval.SumScore.dom2.RESC <- pt(summary(mITEM.dom2.RESC)$coef["TRT.LABLED2.test","t value"],df=ni+nj-3,lower.tail = T)
  
  
  pval.SumScore.dom3.RESC <- pt(summary(mITEM.dom3.RESC)$coef["TRT.LABLED2.test","t value"],df=ni+nj-3,lower.tail = T)
  
  #####
  
  m3.RESC <-  lm(formula=ITEMV1.18 ~ TRT.LABLED + ITEMV1.1, data = dat2.piw.RESC )
  m4.RESC <-  lm(formula=ITEMV2.18 ~ TRT.LABLED + ITEMV2.1, data = dat2.piw.RESC )
  m5.RESC <-  lm(formula=ITEMV3.18 ~ TRT.LABLED + ITEMV3.1, data = dat2.piw.RESC )
  m12.RESC <- lm(formula=ITEMV4.18 ~ TRT.LABLED + ITEMV4.1, data = dat2.piw.RESC )
  m13.RESC <-  lm(formula=ITEMV5.18 ~ TRT.LABLED + ITEMV5.1, data = dat2.piw.RESC )
  m24.RESC <-  lm(formula=ITEMV6.18 ~ TRT.LABLED + ITEMV6.1, data = dat2.piw.RESC )
  m25.RESC <-  lm(formula=ITEMV7.18 ~ TRT.LABLED + ITEMV7.1, data = dat2.piw.RESC )
  m26.RESC <- lm(formula=ITEMV8.18 ~ TRT.LABLED + ITEMV8.1, data = dat2.piw.RESC )
  m27.RESC <-  lm(formula=ITEMV9.18 ~ TRT.LABLED + ITEMV9.1, data = dat2.piw.RESC )
  m28.RESC <-  lm(formula=ITEMV10.18 ~ TRT.LABLED + ITEMV10.1, data = dat2.piw.RESC )
  
  
  MMMall.RESC <- mmm(m3.RESC,m4.RESC,m5.RESC,m12.RESC,m13.RESC,m24.RESC,m25.RESC,m26.RESC,m27.RESC,m28.RESC)
  
  
  contfun.RESC <- function(xx) c(rep(0,3*xx-2),1,rep(0,30-(3*xx-2+1)))
  
  mcontf.RESC <- matrix( c(contfun.RESC(1),contfun.RESC(2),contfun.RESC(3),contfun.RESC(4),contfun.RESC(5),
                      contfun.RESC(6),contfun.RESC(7),contfun.RESC(8),contfun.RESC(9),contfun.RESC(10)), nrow = 10, ncol = 30 ,byrow = TRUE )
  
  
  glhij.RESC <- glht(MMMall.RESC,linfct = mcontf.RESC,alternative="less")
  
  #####
  
  t.eachij.RESC <- coef(glhij.RESC)/sqrt(diag(vcov(glhij.RESC)))
  gencovij.pool.RESC <- vcov(glhij.RESC)
  gencorij.pool.RESC <- cov2cor(gencovij.pool.RESC)
  
  
  pval.eachij.RESC <- numeric(nitem)
  for(k in 1:nitem){
    pval.eachij.RESC[k] <- pt(t.eachij.RESC[k],df=ni+nj-3,lower.tail = T)
  }
  
  
  t.OLSij.RESC <- sum(t.eachij.RESC)/sqrt(sum(gencorij.pool.RESC))
  dfOGij <- 0.5*(ni+nj-3)*(1+1/nitem^2)
  pval.OLSij.RESC <- pt(t.OLSij.RESC,df=dfOGij,lower.tail = T)
  
  
  WOLSij.RESC <- 1/sqrt(sum(gencorij.pool.RESC))
  
  
  WGLSij.RESC <- colSums(solve(gencorij.pool.RESC))/sqrt(sum(solve(gencorij.pool.RESC)))
  names(WGLSij.RESC) <- c(3,4,5,12,13,24,25,26,27,28)
  
  #####
  
  if(length(WGLSij.RESC[WGLSij.RESC<0])>=1){
    negItEM.RESC <- as.integer(names(WGLSij.RESC)[WGLSij.RESC<0])
    NegRij.RESC <- 1
  }else{
    negItEM.RESC <- 100
    NegRij.RESC <- 0
  }
  
  
  t.GLSij.RESC <- sum(solve(gencorij.pool.RESC)%*%t.eachij.RESC)/sqrt(sum(solve(gencorij.pool.RESC)))
  pval.GLSij.RESC <- pt(t.GLSij.RESC,df=dfOGij,lower.tail = T)
  
  #####
  
  gencovij.pool26.RESC <- gencovij.pool.RESC[-c(8),-c(8)]
  gencorij.pool26.RESC <- cov2cor(gencovij.pool26.RESC)
  
  WGLSij26.RESC <- colSums(solve(gencorij.pool26.RESC))/sqrt(sum(solve(gencorij.pool26.RESC)))
  names(WGLSij26.RESC) <- c(3,4,5,12,13,24,25,27,28)
  
  if(length(WGLSij26.RESC[WGLSij26.RESC<0])>=1){
    negItEM26.RESC <- as.integer(names(WGLSij26.RESC)[WGLSij26.RESC<0]) 
    NegRij26.RESC <- 1
  }else{
    negItEM26.RESC <- 100
    NegRij26.RESC <- 0
  }
  
  t.GLSij26.RESC <- sum(solve(gencorij.pool26.RESC)%*%t.eachij.RESC[-8])/sqrt(sum(solve(gencorij.pool26.RESC)))
  
  
  dfOGij26 <- 0.5*(ni+nj-3)*(1+1/(nitem-1)^2)
  pval.GLSij26.RESC <- pt(t.GLSij26.RESC,df=dfOGij26,lower.tail = T)
  res.GLSij26.RESC <- ifelse(pval.GLSij26.RESC<alphaa,1,0)
  
  #####
  
  
  res.OLSij.RESC <- ifelse(pval.OLSij.RESC<alphaa,1,0)
  res.GLSij.RESC <- ifelse(pval.GLSij.RESC<alphaa,1,0)
  
  #####
  
  pval.bonfij.RESC <- nitem*min(pval.eachij.RESC)
  if(pval.bonfij.RESC>=1)  pval.bonfij.RESC <- 1
  res.bonfij.RESC <- ifelse(pval.bonfij.RESC<alphaa,1,0)
  
  #####
  
  pval.Simesij.RESC <- min(hommel::p.adjust(hommel(pval.eachij.RESC)))
  res.Simesij.RESC <- ifelse(pval.Simesij.RESC<alphaa,1,0)
  
  #####
  
  z.min.RESC <- min(qnorm(pt(t.eachij.RESC,df=ni+nj-3)))
  pval.tmin.RESC <- 1-pmvnorm(lower=rep(z.min.RESC,nitem),upper=rep(Inf,nitem),corr = as.matrix(gencorij.pool.RESC) )[1]
  res.tminij.RESC <- ifelse(pval.tmin.RESC<alphaa,1,0)
  
  #####
  
  pval.SumScoreij.RESC <- pt(summary(mITEMall.RESC)$coef["TRT.LABLED2.test","t value"],df=ni+nj-3,lower.tail = T)
  res.SumScoreij.RESC <- ifelse(pval.SumScoreij.RESC<alphaa,1,0)
  
  #####
  
  pval.PSIfij.RESC <- pt(summary(mPSIf.RESC)$coef["TRT.LABLED2.test","t value"],df=ni+nj-3,lower.tail = T)
  res.PSIfij.RESC <- ifelse(pval.PSIfij.RESC<alphaa,1,0)
  
  
  pval.linmPSIbpfij.RESC <- pt(summary(m.linmPSIbpf.RESC)$coef["TRT.LABLED2.test","t value"],df=ni+nj-3,lower.tail = T)
  res.linmPSIbpfij.RESC <- ifelse(pval.linmPSIbpfij.RESC<alphaa,1,0)
  
  #####
  
  pval.Omnibus.RESC <- omnibus.test(pval.eachij.RESC,csStat.H0)
  res.Omnibusij.RESC <- ifelse(pval.Omnibus.RESC<alphaa,1,0)
  
  #####
  
  
  pval.Omnibus3Domij.RESC <- omnibus.test(c(pval.SumScore.dom1.RESC,pval.SumScore.dom2.RESC,pval.SumScore.dom3.RESC),csStat.H03dom)
  res.Omnibus3Domij.RESC <- ifelse(pval.Omnibus3Domij.RESC<alphaa,1,0)
  
  ###########################################
  
  return(list(pval.OLSij=pval.OLSij,res.OLSij=res.OLSij,pval.GLSij=pval.GLSij,res.GLSij=res.GLSij,
              pval.bonfij=pval.bonfij,pval.tmin=pval.tmin,pval.SumScoreij=pval.SumScoreij,res.bonfij=res.bonfij,
              res.tminij=res.tminij,res.SumScoreij=res.SumScoreij,pval.Simesij=pval.Simesij,res.Simesij=res.Simesij,
              pval.Omnibus=pval.Omnibus,res.Omnibusij=res.Omnibusij,res.GLSij26=res.GLSij26,df.residual=df.residual,
              pval.PSIfij=pval.PSIfij,res.PSIfij=res.PSIfij,
              pval.linmPSIbpfij=pval.linmPSIbpfij,res.linmPSIbpfij=res.linmPSIbpfij,
              pval.Omnibus3Domij=pval.Omnibus3Domij,res.Omnibus3Domij=res.Omnibus3Domij,
              pval.OLSij.RESC=pval.OLSij.RESC,res.OLSij.RESC=res.OLSij.RESC,pval.GLSij.RESC=pval.GLSij.RESC,res.GLSij.RESC=res.GLSij.RESC,
              pval.bonfij.RESC=pval.bonfij.RESC,pval.tmin.RESC=pval.tmin.RESC,pval.SumScoreij.RESC=pval.SumScoreij.RESC,res.bonfij.RESC=res.bonfij.RESC,
              res.tminij.RESC=res.tminij.RESC,res.SumScoreij.RESC=res.SumScoreij.RESC,pval.Simesij.RESC=pval.Simesij.RESC,res.Simesij.RESC=res.Simesij.RESC,
              pval.Omnibus.RESC=pval.Omnibus.RESC,res.Omnibusij.RESC=res.Omnibusij.RESC,res.GLSij26.RESC=res.GLSij26.RESC,df.residual.RESC=df.residual.RESC,
              pval.PSIfij.RESC=pval.PSIfij.RESC,res.PSIfij.RESC=res.PSIfij.RESC,
              pval.linmPSIbpfij.RESC=pval.linmPSIbpfij.RESC,res.linmPSIbpfij.RESC=res.linmPSIbpfij.RESC,
              pval.Omnibus3Domij.RESC=pval.Omnibus3Domij.RESC,res.Omnibus3Domij.RESC=res.Omnibus3Domij.RESC))
  
}


#############################################################
# Parallel calculation of simulated power 
#############################################################


### Power simulations using "future_map" function for parallel calculation


###  Adjust the number of cores based on your machine's capabilities to perform parallel calculation
n_cores <- detectCores()


#plan(multicore, workers = n_cores)
plan(multisession, workers = n_cores-1) #server


EFFECTRatio <- seq(0.45,0.75,by=0.05)


Dvec.all.list <- as.list(EFFECTRatio)
names(Dvec.all.list) <- paste0("delta",1:length(EFFECTRatio))

###

start <- Sys.time()

set.seed(123456789)
POWsimSS_list <- future_map(Dvec.all.list, function(effectRatio) {
  
  # Inner loop: Simulation runs for a specific effectRatio
  
  future_map(1:n.sim2, function(i) {
    
    TEST.IRIfun(ni = nmean, nj = nmean, dat = datAbbVie, effectRatio = effectRatio,
                csStat.H0=csStat.H0,csStat.H03dom=csStat.H03dom,
                IRTmodelAbb.pullC = IRTmodelAbb.pullC,
                Linmodelpull.CPSIf = Linmodelpull.CPSIf,
                Linmodelpull.CPSIpf = Linmodelpull.CPSIpf,
                IRTmodelAbb.pullC.RESC = IRTmodelAbb.pullC.RESC,
                Linmodelpull.CPSIf.RESC = Linmodelpull.CPSIf.RESC,
                Linmodelpull.CPSIpf.RESC = Linmodelpull.CPSIpf.RESC)
  }, .options = furrr_options(seed = TRUE))
}, .options = furrr_options(seed = TRUE))



print( Sys.time() - start )


######## Extract the results of the power simulations for the tests based on the simulated data with the original scores

test.names <- c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS", "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom"  ) 
names(POWsimSS_list) <- names(Dvec.all.list)

# Create a list to store the results for each test
results_list_ORIG <- lapply(POWsimSS_list, function(results_delta) {
  
  POWsimSS.rPSIf <- sapply(results_delta, function(x) x$res.PSIfij)
  POWsimSS.rlinmPSIbpf <- sapply(results_delta, function(x) x$res.linmPSIbpfij)
  POWsimSS.rSUMS <- sapply(results_delta, function(x) x$res.SumScoreij)
  POWsimSS.rO <- sapply(results_delta, function(x) x$res.OLSij)
  POWsimSS.rG <- sapply(results_delta, function(x) x$res.GLSij)
  POWsimSS.rG26 <- sapply(results_delta, function(x) x$res.GLSij26)
  POWsimSS.rB <- sapply(results_delta, function(x) x$res.bonfij)
  POWsimSS.rTM <- sapply(results_delta, function(x) x$res.tminij)
  POWsimSS.rSimes <- sapply(results_delta, function(x) x$res.Simesij)
  POWsimSS.rOmnibus <- sapply(results_delta, function(x) x$res.Omnibusij)
  POWsimSS.rOmnibus3Dom <- sapply(results_delta, function(x) x$res.Omnibus3Domij)
  colMeans(cbind(POWsimSS.rPSIf,POWsimSS.rlinmPSIbpf,POWsimSS.rSUMS,POWsimSS.rO,POWsimSS.rG,POWsimSS.rG26,
                 POWsimSS.rB,POWsimSS.rTM,POWsimSS.rSimes,POWsimSS.rOmnibus,POWsimSS.rOmnibus3Dom))
})


# Combine the results into a data frame or matrix
final.res.ORIG <- do.call(rbind, results_list_ORIG)

# Optionally, add column names to the data frame
row.names(final.res.ORIG) <- names(Dvec.all.list)
colnames(final.res.ORIG)<- test.names

final.res.ORIG


######## Extract the results of the power simulations for the tests based on the rescored simulated data

results_list_RESC <- lapply(POWsimSS_list, function(results_delta) {
  
  POWsimSS.rPSIf.RESC <- sapply(results_delta, function(x) x$res.PSIfij.RESC)
  POWsimSS.rlinmPSIbpf.RESC <- sapply(results_delta, function(x) x$res.linmPSIbpfij.RESC)
  POWsimSS.rSUMS.RESC <- sapply(results_delta, function(x) x$res.SumScoreij.RESC)
  POWsimSS.rO.RESC <- sapply(results_delta, function(x) x$res.OLSij.RESC)
  POWsimSS.rG.RESC <- sapply(results_delta, function(x) x$res.GLSij.RESC)
  POWsimSS.rG26.RESC <- sapply(results_delta, function(x) x$res.GLSij26.RESC)
  POWsimSS.rB.RESC <- sapply(results_delta, function(x) x$res.bonfij.RESC)
  POWsimSS.rTM.RESC <- sapply(results_delta, function(x) x$res.tminij.RESC)
  POWsimSS.rSimes.RESC <- sapply(results_delta, function(x) x$res.Simesij.RESC)
  POWsimSS.rOmnibus.RESC <- sapply(results_delta, function(x) x$res.Omnibusij.RESC)
  POWsimSS.rOmnibus3Dom.RESC <- sapply(results_delta, function(x) x$res.Omnibus3Domij.RESC)
  colMeans(cbind(POWsimSS.rPSIf.RESC, POWsimSS.rlinmPSIbpf.RESC,
                 POWsimSS.rSUMS.RESC,POWsimSS.rO.RESC,POWsimSS.rG.RESC,POWsimSS.rG26.RESC,POWsimSS.rB.RESC,POWsimSS.rTM.RESC,
                 POWsimSS.rSimes.RESC,POWsimSS.rOmnibus.RESC,POWsimSS.rOmnibus3Dom.RESC))
})


# Combine the results into a data frame or matrix
final.res.RESC<- do.call(rbind, results_list_RESC)

# Optionally, add column names to the data frame
row.names(final.res.RESC) <- names(Dvec.all.list)
colnames(final.res.RESC)<- test.names

final.res.RESC


#############################################################
# Plot of simulated power of all scenarios and tests 
#############################################################

### Creating the wide format of the results of Original values

POWsim1.ORIG <- data.frame(final.res.ORIG)


POWsim1.ORIG <- POWsim1.ORIG %>%
  mutate(Effect=EFFECTRatio)


POWsim1.ORIG.long <- POWsim1.ORIG %>% 
  pivot_longer(-Effect, names_to = "Test", values_to = "Power")

POWsim1.ORIG.long$Test <- factor(POWsim1.ORIG.long$Test, levels = unique(POWsim1.ORIG.long$Test))


###################################

### Creating the wide format of the results of rescored values

POWsim1.RESC <- data.frame(final.res.RESC)


POWsim1.RESC <- POWsim1.RESC %>%
  mutate(Effect=EFFECTRatio)


POWsim1.RESC.long <- POWsim1.RESC %>% 
  pivot_longer(-Effect, names_to = "Test", values_to = "Power")

POWsim1.RESC.long$Test <- factor(POWsim1.RESC.long$Test, levels = unique(POWsim1.RESC.long$Test))


###################################

POWsim1.ORIG.long$Score <- "Original"
POWsim1.RESC.long$Score <- "Rescore"

POWsim1.BOTH.long <- rbind(POWsim1.ORIG.long,POWsim1.RESC.long)

POWsim1.BOTH.long$Test <- factor(POWsim1.BOTH.long$Test, levels = unique(POWsim1.BOTH.long$Test))


###################################


# Define a function to create legend labels with subscripted values and EFFECTRatio values
create_legend_labels <- function(values, ratios) {
  lapply(seq_along(values), function(i) {
    formatted_ratio <- sprintf("%.2f", ratios[i])
    bquote(rho[.(i)] == .(formatted_ratio))
    
  })
}


legend_labels <- create_legend_labels(1:length(EFFECTRatio), EFFECTRatio)

###################################

y_breaks <- seq(0, 1, by = 0.1)


ggplot(POWsim1.BOTH.long, aes(x = Test, y = Power, color = as.factor(Effect), linetype = Score, group = paste(Effect, Score))) +
  geom_line() +
  geom_point() +
  labs(x = "Test", y = "Power") +
  theme_minimal() +
  scale_color_discrete(name = "Effect", labels = legend_labels) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_y_continuous(breaks = y_breaks, expand = expansion(add = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.title = element_text(size= 14),
        legend.text = element_text(size=11))


