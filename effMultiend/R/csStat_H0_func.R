#' Cumulative mean values of independently generated p-values under the null hypothesis.
#'
#' @description This function generates the independent p-values under the null
#' from the uniform distribution. It then calculates cumulative means of
#' transformed sorted p-values to be used in the implementation of the omnibus test.
#' The transformation as \eqn{h(p)=1/p} was considered to give weight to extreme p-values.
#'
#'
#' @param m Number of individual test statistics (hypothesis tests) in the global test.
#' @param N.sim Number of simulations for calculation of cumulative mean
#' values of generated p-values.
#'
#' @return A matrix with the dimension `m*N.sim` containing the cumulative mean values
#' of individual p-values which are possibly transformed and then sorted.
#' Each column corresponds to a simulation.
#'
#'
#' @seealso \code{\link{omnibus_test}}
#'
#' @importFrom stats runif
#' @importFrom dplyr cummean
#' @importFrom memoise memoise
#'
#' @references
#' Futschik, A., Taus, T., & Zehetmayer, S. (2019). An omnibus test for the global
#' null hypothesis. \emph{Statistical methods in medical research}, **28**(8), 2292-2304.
#'
#' @author Sonja Zehetmayer
#'

csStat_H0_func<- memoise::memoise(function(m,N.sim){

  #if(!is.null(seed)) set.seed(seed)
  stat.H0.p <- matrix(runif(N.sim*m), nrow = m)
  stat.H0 <- 1/stat.H0.p # the data is (reciprocal) transformed
  stat.H0sort <- apply(stat.H0,2,sort,decreasing=TRUE) #sorting
  csStat.H0 <- apply(stat.H0sort,2, dplyr::cummean) # cumulative mean values are calculated
  csStat.H0
})

