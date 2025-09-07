#' The Omnibus test for testing the global null hypothesis based on independent p-values
#'
#' @description
#' The function calculates an Omnibus test statistic using the given independent
#' p-values and the null distribution of cumulative mean values. It focuses on
#' a general situation where independent p-values are available from several hypothesis tests that
#' are assumed to be uniformly distributed under the null hypothesis. It is an
#' alternative test for the global null hypothesis of no treatment effect in any
#' item which is based on cumulative sums of (possibly) transformed sorted p-values.
#'
#' @param p A vector of p-values.
#' @param csStat.H00 Null distribution of cumulative mean values.
#'
#' @return The omnibus test statistic.
#'
#' @seealso \code{\link{csStat_H0_func}}
#'
#'
#' @importFrom stats ecdf
#' @importFrom matrixStats colCumsums
#' @importFrom matrixStats colMaxs
#'
#'
#' @references
#' Futschik, A., Taus, T., & Zehetmayer, S. (2019). An omnibus test for the
#' global null hypothesis. \emph{Statistical methods in medical research}, **28**(8), 2292-2304.
#'
#' @author Sonja Zehetmayer
#'



omnibus_test  <- function (p, csStat.H00)
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
