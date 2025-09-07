#' Estimation of an Item Response Theory (IRT) model, a linear model of
#' items and the corresponding latent traits.
#'
#' @description
#' This function receives a simulated data set in the long format as the input.
#' The input simulated data set resembles items scores of the Progressive
#' Supranuclear Palsy (PSP) as in the ABBV-8E12 trial. The simulated data set
#' must contain specific variables as explained in data description part, below.
#' After accumulating data from independent time points (visits, \code{VISIT})
#' with the aim to enhance the precision of the model estimates, the graded
#' response IRT model is estimated. Using the estimated latent trait (progression rate)
#' for each patient from the IRT model, as the dependent variable, the function
#' fits a linear model with the 10 FDA recommended items as independent variables.
#' This yields an approximation of the latent variable based on a weighted sum
#' of the item scores. This novel approach for approximation of the latent trait
#' is simpler to compute without going through complex derivations. Similar
#' estimations are performed after the items are rescored based on the FDA suggestion.
#' In addition to the above estimations, the function also returns the patient IDs
#' who were present at both time points of interest (baseline and week 52
#' corresponding to \code{VISIT=1}, \code{VISIT=18}).
#'
#' @param longdata A data frame containing repeated measures of patients.
#' This dataset is used for estimation of the IRT model, linear model of items,
#' and corresponding latent trait estimation.
#'   The variables include:
#'     - \item{ID}{Unique patient identifier}
#'     - \item{STUDY_ID}{1, referring to the ABBV-8E12 trial}
#'     - \item{AGE}{Age of patients}
#'     - \item{Sex}{1 for Female, 2 for Male}
#'     - \item{PHENOTYPE}{PSP-RS ~ 9, PSP Diagnostic phenotype}
#'     - \item{TRT}{Treatment, simulated treatment labels, 1 corresponds to the test group and 0 to the control group}
#'     - \item{VISIT}{Visit number, where VISIT can take the values 1, 6, 11, 14, 18. VISIT=1 and VISIT=18 correspond to baseline and week 52}
#'     - \item{TIME}{Time in year}
#'     - \item{DV}{Dependent variable, scores of the 10-item version of the PSPRS recommended by FDA, mostly vary from 0 to 4}
#'     - \item{ITEM}{Item measured, representing the 10-item version of the PSPRS recommended by the FDA. Originally, these items were numbered 3, 4, 5, 12, 13, 24, 25, 26, 27, 28. In our implementation, we renamed them from 1 to 10}
#'     - \item{DOMAIN}{Domain number}
#'     - \item{SCALE}{PSPRS ~ 1, Global scale}
#'     - \item{RICH}{Indicator variable, taking the value 1 if PHENOTYPE is 9, otherwise 0}
#'     - \item{INTERV}{1, if the study has been intervention}
#'
#' @return A list containing various model estimates with the following components.
#' Note that values with the extension '.RESC' stand for the rescored items and have similar interpretations as those in the following:
#' \item{IRTmodel.full}{Full IRT model}
#' \item{IRTpar.est}{Estimated parameters for the IRT model.}
#' \item{Linmodel.PSIpf}{Linear model of items with the (transformed) estimated
#' latent trait from the IRT model as the depedent variable}
#' \item{est.IRTlatent.pf}{Estimated latent traits for each patient from the IRT model}
#' \item{est.LMlatent.bpf}{Approximated latent traits for each patient from the linear model of items}
#' \item{id.overlap}{Patient IDs who have all item scores both at baseline and week 52}

#'
#'
#' @details
#'  The input scores include the (simulated) 10-item version of the PSPRS recommended by FDA.
#'  The item scores are categorical values which mostly vary from zero to 4. The estimations based on the
#'  true data set of the ABBV-8E12 trial are presented in Yousefi et al. (2023).
#'
#' @importFrom stats complete.cases plogis lm predict qlogis fitted
#' @importFrom mirt coef
#' @importFrom mirt mirt
#' @importFrom dplyr "%>%"
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_wider
#' @importFrom readr read_csv
#'
#'
#' @references
#' Ueckert, S. (2018). Modeling composite assessment data using item response theory.
#' \emph{CPT: Pharmacometrics & Systems Pharmacology}, **7**(4), 205-218.
#' Yousefi, E., Gewily, M., K\"{o}nig, F., H\"{o}glinger, G., Hopfner, F.,
#' Karlsson, M. O., Ristl, R., Zehetmayer, S. & Posch, M. (2023).
#' Efficiency of Multivariate Tests in Trials in Progressive Supranuclear Palsy.
#' \emph{arXiv preprint arXiv:2312.08169}.
#'
#' @author Elham Yousefi
#'
#' @seealso
#' \code{\link{simPSPdata}}
#'
#' @examples
#' # Example usage:
#' \dontrun{
#' sim_data <- read.csv("C:/2024/PSP_directory/MultiendPSP/effMultiend/inst/extdata/simPSP_lim.csv")
#' result <- est_funcs(sim_data)
#' }
#'
#' @export


#test example
#library(tidyr)
#library(dplyr)
#library(mirt)
#sim_data<- read.csv("C:/2024/PSP_directory/MultiendPSP/effMultiend/inst/extdata/simPSP_lim.csv")
#sim_data<- read.csv("C:/2024/MUWfiles/2023_08_PSP_DD/artificial_data.csv")

est_funcs <- function(longdata){


  #browser()


  # Check if required variables are present in the input data
  required.vars <- c("VISIT", "ITEM", "ID", "TRT", "DV")
  if (any(!required.vars %in% colnames(longdata))) {
  stop("Input data is missing one or more required variables: VISIT, ITEM, ID, TRT, DV")
  }


  if ( !all(longdata$ID %% 1 == 0) || any(longdata$ID < 1) || !is.finite(max(longdata$ID))) {
  stop("Invalid values in the ID column. IDs should be integers starting from one with a finite maximum.")
  }


  if (!all(longdata$TRT %in% c(0, 1))) {
  stop("Invalid values in the TRT column. Treatment labels should be either 0 or 1.")
  }


  if (!all(longdata$ITEM %in% 1:10)) {
  stop("Invalid values in the ITEM column. ITEMs must vary from 1 to 10.")
  }



  if (!all(longdata$DV %in% 0:4)) {
  stop("Invalid values in the DV column. DV must have values from 0 to 4.")
  }


  required_visits <- c(1, 6, 11, 14, 18)
  if (any(!required_visits %in% unique(longdata$VISIT))) {
  stop("Input data is missing one or more required visit numbers: 1, 6, 11, 14, 18")
  }

  ###

  id.summary1 <- longdata %>% group_by(ID) %>% filter(VISIT==18) %>% pull(ID) %>% unique()


  id.summary2 <- longdata %>% group_by(ID) %>% filter(VISIT==1) %>% pull(ID) %>% unique()

  # intersection of the patient IDs
  id.overlap <- intersect(id.summary1,id.summary2)

  ###

  dataVISITS.wide <- longdata %>% filter(VISIT %in% c(1,6,11,14,18)) %>%
    pivot_wider(id_cols = c(ID,TRT), names_from = c(ITEM,VISIT), values_from = DV, names_sep = ".", names_prefix = "ITEMV")



  ### select the data at baseline (only complete cases where the data for all items is available) so that residuals and fitted values can be evaluated in further steps
  dataVISIT1.wide <- dataVISITS.wide %>% dplyr::select(ends_with(".1"))
  dataVISIT1.wideC <- dataVISIT1.wide[complete.cases(dataVISIT1.wide),]

  colnames(dataVISIT1.wideC) <- paste0("ITEM",1:10)



  ### select the data at visit 6 (only complete cases where the data for all items is available), similarly
  dataVISIT6.wide <- dataVISITS.wide %>% dplyr::select(ends_with(".6"))
  dataVISIT6.wideC <- dataVISIT6.wide[complete.cases(dataVISIT6.wide),]

  colnames(dataVISIT6.wideC) <- paste0("ITEM",1:10)



  ### select the data at visit 11 (only complete cases where the data for all items is available), similarly
  dataVISIT11.wide <- dataVISITS.wide %>% dplyr::select(ends_with(".11"))
  dataVISIT11.wideC <- dataVISIT11.wide[complete.cases(dataVISIT11.wide),]

  colnames(dataVISIT11.wideC) <- paste0("ITEM",1:10)



  ### select the data at visit 14 (only complete cases where the data for all items is available), similarly
  dataVISIT14.wide <- dataVISITS.wide %>% dplyr::select(ends_with(".14"))
  dataVISIT14.wideC <- dataVISIT14.wide[complete.cases(dataVISIT14.wide),]

  colnames(dataVISIT14.wideC) <- paste0("ITEM",1:10)



  ### select the data at visit 18 (only complete cases where the data for all items is available), similarly
  dataVISIT18.wide <- dataVISITS.wide %>% dplyr::select(ends_with(".18"))
  dataVISIT18.wideC <- dataVISIT18.wide[complete.cases(dataVISIT18.wide),]

  colnames(dataVISIT18.wideC) <- paste0("ITEM",1:10)



  ### combining/pooling the selected data from visits of interest to enhance the precision of the IRT model estimates
  datapool.wideC <- rbind(dataVISIT1.wideC,dataVISIT6.wideC,dataVISIT11.wideC,dataVISIT14.wideC,dataVISIT18.wideC)



  ### Check the unique values of each item before rescoring
  #apply(datapool.wideC, 2, function(x) unique(x))



  ### Rescoring the items based on the FDA recommendation. Note that items 3 and 13, recoded as 1 and 5 in implementations, don't contain any changes.

  datapool.wideC.RESC <- datapool.wideC %>%
    dplyr::mutate(
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
  #apply(datapool.wideC.RESC, 2, function(x) unique(x))



  ### Fit the graded response IRT model to the dataset of the ABBV-8E12 trial, where aggregated data with original scores from all treatment groups and visits are used for estimation
  catg.item <- as.vector(apply(datapool.wideC, 2, function(x) length(unique(x)))) # get the number of categories for each item

  IRTmodel.full <- mirt(data=datapool.wideC,model=1, itemtype = "graded",GR = catg.item, verbose = 0) # fitting the IRT model

  IRTpar.est <- coef(IRTmodel.full,simplify=T,IRTpars=T) # estimation of the IRT model parameters

  ###

  ### Approximation of the latent variable based on a linear model fit, where as above the aggregated data from ABBV-8E12 trial are considered
  ### The dependent variable, estimated latent variable from the IRT fitted model, is transformed using the standard logistic function to improve model fit
  est.IRTlatent.pf <- plogis(fscores(IRTmodel.full))
  datapool.wideC.PSI <- datapool.wideC %>% mutate(est.IRTlatent.pf=est.IRTlatent.pf)

  Linmodel.PSIpf <- lm(est.IRTlatent.pf ~ ITEM1 + ITEM2 + ITEM3 + ITEM4 +ITEM5 + ITEM6 +ITEM7 +
                         ITEM8 +ITEM9 + ITEM10, data = datapool.wideC.PSI)


  fitted.LM <- predict(Linmodel.PSIpf, newdata = datapool.wideC)
  fitted.LM <- pmin(pmax(0.001,fitted.LM),0.999) # out of bound values are set to 0.001-0.999 (close to 0-1)

  est.LMlatent.bpf <- qlogis(fitted.LM)


  ######### Similar estimation of the graded response IRT model and the linear model fit for the rescored aggregated data from ABBV-8E12 trial

  catg.item.RESC <- as.vector(apply(datapool.wideC.RESC, 2, function(x) length(unique(x)))) # get the number of categories for each item

  IRTmodel.full.RESC <- mirt(data=datapool.wideC.RESC,model=1, itemtype = "graded",GR = catg.item.RESC, verbose = 0) # fitting the IRT model

  IRTpar.est.RESC <- coef(IRTmodel.full.RESC,simplify=T,IRTpars=T) # estimation of the IRT model parameters

  ###

  est.IRTlatent.pf.RESC <- plogis(fscores(IRTmodel.full.RESC))
  datapool.wideC.RESC.PSI <- datapool.wideC.RESC %>% mutate(est.IRTlatent.pf.RESC=est.IRTlatent.pf.RESC)


  Linmodel.PSIpf.RESC <- lm(est.IRTlatent.pf.RESC ~ ITEM1 + ITEM2 + ITEM3 + ITEM4 +ITEM5 + ITEM6 +ITEM7 +
                              ITEM8 +ITEM9 + ITEM10, data = datapool.wideC.RESC.PSI)


  fitted.LM.RESC <- predict(Linmodel.PSIpf.RESC, newdata = datapool.wideC.RESC)
  fitted.LM.RESC <- pmin(pmax(0.001,fitted.LM.RESC),0.999)

  est.LMlatent.bpf.RESC <- qlogis(fitted(Linmodel.PSIpf.RESC))

  ###

  return(list(IRTmodel.full=IRTmodel.full, IRTpar.est=IRTpar.est, Linmodel.PSIpf=Linmodel.PSIpf,id.overlap=id.overlap,
              IRTmodel.full.RESC=IRTmodel.full.RESC, IRTpar.est.RESC=IRTpar.est.RESC, Linmodel.PSIpf.RESC=Linmodel.PSIpf.RESC,
              est.IRTlatent.pf=est.IRTlatent.pf,est.IRTlatent.pf.RESC=est.IRTlatent.pf.RESC,
              est.LMlatent.bpf=est.LMlatent.bpf,est.LMlatent.bpf.RESC=est.LMlatent.bpf.RESC))



}



