#' parallel run of the IRT_based Simulations
#'
#'
#' @description
#' This function runs in parallel the IRT based simulations for estimation of
#' the type one error rate and simulation of the power of statistical hypothesis
#' tests. The analysis methods are described in \code{\link{test_fun_IRT}}.
#' If \code{effectRatio.list=NULL} the simulations are run to estimate the type one error rate.
#'
#' @param ni Number of subjects in the test group.
#' @param nj Number of subjects in the control group.
#' @param dat A data set in the long format, containing the necessary variables
#' as in \code{\link{simPSPdata}} data set to perform the simulations.
#' @param effectRatio.list A list of numeric values or NULL. A list of the rho
#'  values, indicating a beneficial treatment effect, to be used in power
#'  simulations. It is used to model disease progression under the alternative.
#'  If NULL, simulations are performed to estimate the type I error rate,
#'  otherwise the power simulations are performed. Default = NULL.
#' @param alphaa Alpha level. Default is 0.025.
#' @param n.sim Number of simulation runs.
#' @param test_names A character vector specifying the names of tests to be performed. This includes all the following tests:
#'   \itemize{
#'     \item \code{"IRT.PSIF"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"LM.PSIBPF"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"SumS"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"OLS"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"GLS"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"GLS_26"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"Bonf"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"Tmin"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"Simes"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"Omnibus"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'     \item \code{"Omnibus_dom"}: refer to the documentation of the \link{test_fun_IRT} function for more details.
#'   }
#' @param csStat.H0 Cumulative mean values of independently generated item-wise p-values under the null hypothesis.
#'     This variable is pre-computed using \code{\link{csStat_H0_func}}. Item-wise p-values are of length 10.
#' @param csStat.H03dom Cumulative mean values of independently generated domain-wise p-values under the null hypothesis.
#'     This variable is pre-computed using \code{\link{csStat_H0_func}}. Domain-wise p-values are of length 3.
#'
#' @return A data frame containing the results of simulations for all the
#' analysis tests i.e., the proportions of time the H0 has been rejected.
#' This includes the results for the original and rescored simulated samples as the following:
#' \item{final.res.ORIG}{Results of simulations of all analysis tests
#'  for the simulated items with the original scores}
#' \item{final.res.RESC}{Results of simulations of all analysis tests
#'  for the rescored simulated items}
#'
#'
#'
#' @importFrom furrr future_map
#' @importFrom furrr furrr_options
#' @importFrom parallel detectCores
#' @importFrom dplyr "%>%"
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_wider
#' @importFrom future plan
#' @importFrom future multisession
#' @import effMultiend
#'
#' @seealso
#' \code{\link{simPSPdata}}
#' \code{\link{test_fun_IRT}} for details on the underlying analysis methods.
#' \code{\link{csStat_H0_func}}
#' \code{\link{omnibus_test}}
#'
#' @references
#' Yousefi, E., Gewily, M., K\"{o}nig, F., H\"{o}glinger, G., Hopfner, F.,
#' Karlsson, M. O., Ristl, R., Zehetmayer, S. & Posch, M. (2023).
#' Efficiency of Multivariate Tests in Trials in Progressive Supranuclear Palsy.
#'  \emph{arXiv preprint arXiv:2312.08169}.
#'
#'
#' @author Elham Yousefi
#'
#' @examples
#' library(furrr)
#' library(parallel)
#' library(dplyr)
#' library(tidyr)
#' library(future)
#' library(matrixStats)
#' library(mirt)
#' library(dplyr)
#' library(tidyr)
#' library(tidyselect)
#' library(Matrix)
#' library(hommel)
#' library(multcomp)
#' library(mvtnorm)
#' \dontrun{
#' # Example usage:
#' Dvec.list <- as.list(seq(0.45,0.75,by=0.05))
#' sim_data <- read.csv("C:/2024/PSP_directory/MultiendPSP/effMultiend/inst/extdata/simPSP_lim.csv")
#' set.seed(2403)
#' csStat.H0item <- csStat_H0_func(m = 10, N.sim = 100000)
#' set.seed(2310)
#' csStat.H0domain <-csStat_H0_func(m=3,N.sim=100000)
#' analys_sim_result <- sim_run_par(ni = 50, nj = 50, dat = sim_data,
#' effectRatio.list = Dvec.list, n.sim = 10000,
#' test_names=c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS", "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom"),
#' csStat.H0=csStat.H0item, csStat.H03dom=csStat.H0domain)
#' }
#'
#' @export



# test example
#library(furrr)
#library(parallel)
#library(dplyr)
#library(tidyr)
#library(future)
#library(GMCM)
#library(matrixStats)
#library(mirt)
#library(dplyr)
#library(tidyr)
#library(tidyselect)
#library(Matrix)
#library(hommel)
#library(multcomp)
#library(mvtnorm)

#Dvec.list <- as.list(seq(0.45,0.75,by=0.05))
#sim_data <- read.csv("C:/2024/PSP_directory/MultiendPSP/effMultiend/inst/extdata/simPSP_lim.csv")
#set.seed(2403)
#csStat.H0item <- csStat_H0_func(m = 10, N.sim = 100000)
#set.seed(2310)
#csStat.H0domain <-csStat_H0_func(m=3,N.sim=100000)

#analys_sim_result <- sim_run_par(ni = 50, nj = 50, dat = sim_data,
#effectRatio.list = Dvec.list, alphaa = 0.025, usefitted=TRUE, n.sim = 100,
#test_names=c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS",
#             "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom"),
#csStat.H0=csStat.H0item, csStat.H03dom=csStat.H0domain)


sim_run_par <- function(ni, nj, dat, effectRatio.list = NULL, alphaa = 0.025, usefitted=TRUE, n.sim,
                        test_names=c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS",
                       "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom")) {

  if(is.null(effectRatio.list)){


    n_cores <- detectCores()
    plan(multisession, workers = n_cores-1)

    start <- Sys.time()

    asimSS_list <- future_map(1:n.sim, function(i) {

      set.seed(i)  # Set a different seed for each simulation run
      test_fun_IRT(ni = ni, nj = nj, dat = dat, effectRatio = effectRatio.list, alphaa = alphaa, usefitted = usefitted)

    }, .options = furrr_options(seed = TRUE))


    print( Sys.time() - start )


    ### Extract the results of the type ore error rate, from the list

    alphasimSS.rO <- mean(sapply(asimSS_list, function(x) x$res.OLSij))
    alphasimSS.rG <- mean(sapply(asimSS_list, function(x) x$res.GLSij))
    alphasimSS.rB <- mean(sapply(asimSS_list, function(x) x$res.bonfij))
    alphasimSS.rTM <- mean(sapply(asimSS_list, function(x) x$res.tmaxij))
    alphasimSS.rSUMS <- mean(sapply(asimSS_list, function(x) x$res.SumScoreij))
    alphasimSS.rPSIf <- mean(sapply(asimSS_list, function(x) x$res.PSIfij))
    alphasimSS.rlinmPSIbpf <- mean(sapply(asimSS_list, function(x) x$res.linmPSIbpfij))
    alphasimSS.rG26 <- mean(sapply(asimSS_list, function(x) x$res.GLSij26))
    alphasimSS.rSimes <- mean(sapply(asimSS_list, function(x) x$res.Simesij))
    alphasimSS.rOmnibus <- mean(sapply(asimSS_list, function(x) x$res.Omnibusij))
    alphasimSS.rOmnibus3Dom <- mean(sapply(asimSS_list, function(x) x$res.Omnibus3Domij))


    #test_names <- c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS", "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom"  )
    test_names<- test_names

    final.res.ORIG <- c(alphasimSS.rPSIf,alphasimSS.rlinmPSIbpf,alphasimSS.rSUMS,alphasimSS.rO,alphasimSS.rG,alphasimSS.rG26,
                        alphasimSS.rB,alphasimSS.rTM,alphasimSS.rSimes,alphasimSS.rOmnibus, alphasimSS.rOmnibus3Dom)

    names(final.res.ORIG) <- test_names

    ###

    alphasimSS.rO.RESC <- mean(sapply(asimSS_list, function(x) x$res.OLSij.RESC))
    alphasimSS.rG.RESC <- mean(sapply(asimSS_list, function(x) x$res.GLSij.RESC))
    alphasimSS.rB.RESC <- mean(sapply(asimSS_list, function(x) x$res.bonfij.RESC))
    alphasimSS.rTM.RESC <- mean(sapply(asimSS_list, function(x) x$res.tmaxij.RESC))
    alphasimSS.rSUMS.RESC <- mean(sapply(asimSS_list, function(x) x$res.SumScoreij.RESC))
    alphasimSS.rPSIf.RESC <- mean(sapply(asimSS_list, function(x) x$res.PSIfij.RESC))
    alphasimSS.rlinmPSIbpf.RESC <- mean(sapply(asimSS_list, function(x) x$res.linmPSIbpfij.RESC))
    alphasimSS.rG26.RESC <- mean(sapply(asimSS_list, function(x) x$res.GLSij26.RESC))
    alphasimSS.rSimes.RESC <- mean(sapply(asimSS_list, function(x) x$res.Simesij.RESC))
    alphasimSS.rOmnibus.RESC <- mean(sapply(asimSS_list, function(x) x$res.Omnibusij.RESC))
    alphasimSS.rOmnibus3Dom.RESC <- mean(sapply(asimSS_list, function(x) x$res.Omnibus3Domij.RESC))


    final.res.RESC <- c(alphasimSS.rPSIf.RESC,alphasimSS.rlinmPSIbpf.RESC,alphasimSS.rSUMS.RESC,alphasimSS.rO.RESC,alphasimSS.rG.RESC,alphasimSS.rG26.RESC,
                        alphasimSS.rB.RESC,alphasimSS.rTM.RESC,alphasimSS.rSimes.RESC,alphasimSS.rOmnibus.RESC, alphasimSS.rOmnibus3Dom.RESC)

    names(final.res.RESC) <- test_names


  }else{

    n_cores <- detectCores()
    plan(multisession, workers = n_cores - 1)

    start <- Sys.time()

    #set.seed(123456789)
    POWsimSS_list <- future_map(effectRatio.list, function(effectRatio) {

      # Inner loop: Simulation runs for a specific effect ratio
      future_map(1:n.sim, function(i) {

        test_fun_IRT(ni = ni, nj = nj, dat = dat, effectRatio = effectRatio, alphaa = alphaa, usefitted = usefitted)
      }, .options = furrr_options(seed = TRUE))

    }, .options = furrr_options(seed = TRUE))

    print( Sys.time() - start )

    test_names <- c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS", "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom")

    names(POWsimSS_list) <- paste0("rho", 1:length(effectRatio.list))


    results_list_ORIG <- lapply(POWsimSS_list, function(results_rho) {

      POWsimSS.rPSIf <- sapply(results_rho, function(x) x$res.PSIfij)
      POWsimSS.rlinmPSIbpf <- sapply(results_rho, function(x) x$res.linmPSIbpfij)
      POWsimSS.rSUMS <- sapply(results_rho, function(x) x$res.SumScoreij)
      POWsimSS.rO <- sapply(results_rho, function(x) x$res.OLSij)
      POWsimSS.rG <- sapply(results_rho, function(x) x$res.GLSij)
      POWsimSS.rG26 <- sapply(results_rho, function(x) x$res.GLSij26)
      POWsimSS.rB <- sapply(results_rho, function(x) x$res.bonfij)
      POWsimSS.rTM <- sapply(results_rho, function(x) x$res.tmaxij)
      POWsimSS.rSimes <- sapply(results_rho, function(x) x$res.Simesij)
      POWsimSS.rOmnibus <- sapply(results_rho, function(x) x$res.Omnibusij)
      POWsimSS.rOmnibus3Dom <- sapply(results_rho, function(x) x$res.Omnibus3Domij)
      colMeans(cbind(POWsimSS.rPSIf,POWsimSS.rlinmPSIbpf,POWsimSS.rSUMS,POWsimSS.rO,POWsimSS.rG,POWsimSS.rG26,
                     POWsimSS.rB,POWsimSS.rTM,POWsimSS.rSimes,POWsimSS.rOmnibus,POWsimSS.rOmnibus3Dom))
    })


    # Combine the results into a data frame or matrix
    final.res.ORIG <- do.call(rbind, results_list_ORIG)

    row.names(final.res.ORIG) <- names(POWsimSS_list)
    colnames(final.res.ORIG)<- test_names


    ### Extract the results of the power simulations for the tests based on the rescored simulated data

    results_list_RESC <- lapply(POWsimSS_list, function(results_rho) {

      POWsimSS.rPSIf.RESC <- sapply(results_rho, function(x) x$res.PSIfij.RESC)
      POWsimSS.rlinmPSIbpf.RESC <- sapply(results_rho, function(x) x$res.linmPSIbpfij.RESC)
      POWsimSS.rSUMS.RESC <- sapply(results_rho, function(x) x$res.SumScoreij.RESC)
      POWsimSS.rO.RESC <- sapply(results_rho, function(x) x$res.OLSij.RESC)
      POWsimSS.rG.RESC <- sapply(results_rho, function(x) x$res.GLSij.RESC)
      POWsimSS.rG26.RESC <- sapply(results_rho, function(x) x$res.GLSij26.RESC)
      POWsimSS.rB.RESC <- sapply(results_rho, function(x) x$res.bonfij.RESC)
      POWsimSS.rTM.RESC <- sapply(results_rho, function(x) x$res.tmaxij.RESC)
      POWsimSS.rSimes.RESC <- sapply(results_rho, function(x) x$res.Simesij.RESC)
      POWsimSS.rOmnibus.RESC <- sapply(results_rho, function(x) x$res.Omnibusij.RESC)
      POWsimSS.rOmnibus3Dom.RESC <- sapply(results_rho, function(x) x$res.Omnibus3Domij.RESC)
      colMeans(cbind(POWsimSS.rPSIf.RESC, POWsimSS.rlinmPSIbpf.RESC,
                     POWsimSS.rSUMS.RESC,POWsimSS.rO.RESC,POWsimSS.rG.RESC,POWsimSS.rG26.RESC,POWsimSS.rB.RESC,POWsimSS.rTM.RESC,
                     POWsimSS.rSimes.RESC,POWsimSS.rOmnibus.RESC,POWsimSS.rOmnibus3Dom.RESC))
    })


    # Combine the results into a data frame or matrix
    final.res.RESC<- do.call(rbind, results_list_RESC)

    row.names(final.res.RESC) <- names(POWsimSS_list)
    colnames(final.res.RESC)<- test_names

  }


  return(list(final.res.ORIG = final.res.ORIG, final.res.RESC = final.res.RESC))
}



