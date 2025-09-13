#' Plot of results of analysis tests of the IRT-based simulations
#'
#' @description
#' This function plots the results of multiple analysis tests derived from the
#' parallel simulations described in \code{\link{sim_run_par}} and
#' \code{\link{test_fun_IRT}}. If \code{effectRatio.=NULL} the estimated
#' type one error rate is plotted. Otherwise, based on the values of
#' \code{effectRatio=unlist(effectRatio.list)} result of the power simulations
#' of analysis tests is plotted.
#'
#'
#' @param simresult Result of simulations for both original and rescored
#' simulated items. See \code{\link{sim_run_par}} for more details. This include the following items:
#'   \itemize{
#'     \item \code{final.res.ORIG}: refer to returned values of the \link{sim_run_par} function for more details. This part of results can be accessed using the dollar operator as simresult$final.res.ORIG.
#'     \item \code{final.res.RESC}: refer to returned values of the \link{sim_run_par} function for more details. This part of results can be accessed using the dollar operator as simresult$final.res.RESC.
#'     }
#' @param effectRatio Vector of the rho values, indicating a beneficial
#' treatment effect, to be used in power simulations. It is used to model disease
#' progression under the alternative. See section 4.1 of the following paper.
#' @param alphaa Alpha level (default is 0.025).
#' @param test.names refer to the \link{sim_run_par} function for more details.
#' @param scoretype A character vector specifying type of scores. \code{scoretype="Original"} refers to the results when the original scaling provided by the FDA has been used.
#' \code{scoretype="Rescore"} refers to the case with the FDA rescoring. See the referenced paper below for more details.
#'
#' @return Plot of the simulation results.
#'
#' @importFrom ggplot2 aes geom_line geom_point labs theme_minimal expansion
#' @importFrom ggplot2 scale_color_discrete scale_y_continuous geom_hline
#' @importFrom ggplot2 coord_cartesian theme element_text scale_linetype_manual
#' @importFrom dplyr "%>%"
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#' @import effMultiend
#'
#' @seealso
#' \code{\link{simPSPdata}},
#' \code{\link{test_fun_IRT}} for details on the underlying analysis methods and
#' \code{\link{sim_run_par}}.
#'
#' @references
#' Yousefi, E., Gewily, M., K\"{o}nig, F., H\"{o}glinger, G., Hopfner, F.,
#' Karlsson, M. O., Ristl, R., Zehetmayer, S. & Posch, M. (2023).
#' Efficiency of Multivariate Tests in Trials in Progressive Supranuclear Palsy.
#' \emph{arXiv preprint arXiv:2312.08169}.
#'
#'
#' @author Elham Yousefi
#'
#' @examples
#' # Example usage:
#' \dontrun{
#' # Example usage:
#' \dontrun{
#' # Load built-in example dataset
#' data("simPSPdata")
#'
#' # Define effect ratios to simulate treatment effects
#' Dvec <- seq(0.45, 0.75, by = 0.05)
#'
#' # Run parallel simulations
#' simresult <- sim_run_par(
#'   ni = 50, nj = 50,
#'   dat = simPSPdata,
#'   effectRatio.list = as.list(Dvec),
#'   n.sim = 1000,  # reduced for example speed
#'   test_names = c("IRT.PSIF", "LM.PSIBPF", "SumS", "OLS",
#'                  "GLS", "GLS_26", "Bonf", "Tmin",
#'                  "Simes", "Omnibus", "Omnibus_dom")
#' )
#'
#' # Plot simulation results
#' plot_IRTsim(
#'   simresult = simresult,
#'   effectRatio = Dvec,
#'   alphaa = 0.025,
#'   test.names = c("IRT.PSIF", "LM.PSIBPF", "SumS", "OLS",
#'                  "GLS", "GLS_26", "Bonf", "Tmin",
#'                  "Simes", "Omnibus", "Omnibus_dom"),
#'   scoretype = c("Original", "Rescore")
#' )
#' }
#'
#'
#' @export



plot_IRTsim <- function(simresult, effectRatio = NULL, alphaa = 0.025,
                        test.names=c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS", "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom"),
                        scoretype=c("Original","Rescore")){

  if(is.null(effectRatio)){


    #test.names <- c("IRT.PSIF", "LM.PSIBPF", "SumS",  "OLS", "GLS", "GLS_26", "Bonf", "Tmin", "Simes", "Omnibus", "Omnibus_dom")


    Alphasim1.ORIG <- data.frame(alphaa.sim=simresult$final.res.ORIG, Test=test.names)
    Alphasim1.RESC <- data.frame(alphaa.sim=simresult$final.res.RESC, Test=test.names)

    rou.up.max <- round((max(Alphasim1.ORIG$alphaa.sim)+0.005),digits = 2)


    y_breaks <- seq(0, rou.up.max, by = 0.005)

    Alphasim1.ORIG$Score <- scoretype[1]
    Alphasim1.RESC$Score <- scoretype[2]

    Alphasim1.BOTH.long <- rbind(Alphasim1.ORIG,Alphasim1.RESC)

    Alphasim1.BOTH.long$Test <- factor(Alphasim1.BOTH.long$Test, levels = unique(Alphasim1.BOTH.long$Test))

    #####################################


    ggplot(Alphasim1.BOTH.long, aes(x = Test, y = alphaa.sim, color = as.factor(Score), group = paste(Score))) +
      geom_line() +
      geom_point() +
      labs(x = "Test", y = "Type I error") +
      theme_minimal() +
      scale_color_discrete(name = "Score") +
      scale_y_continuous(breaks = y_breaks) +
      geom_hline(yintercept = alphaa, linetype = "dashed") +
      coord_cartesian(ylim = c(0,rou.up.max)) +
      theme(axis.text = element_text(size = 11),
            axis.title = element_text(size = 14),
            legend.title = element_text(size= 14),
            legend.text = element_text(size=11))



  }else{

    EFFECTRatio <- effectRatio
    final.res.orig <- simresult$final.res.ORIG

    POWsim1.ORIG <- data.frame(final.res.orig)


    POWsim1.ORIG <- POWsim1.ORIG %>%
      mutate(Effect=EFFECTRatio)


    POWsim1.ORIG.long <- POWsim1.ORIG %>%
      pivot_longer(-Effect, names_to = "Test", values_to = "Power")

    POWsim1.ORIG.long$Test <- factor(POWsim1.ORIG.long$Test, levels = unique(POWsim1.ORIG.long$Test))

    ##
    final.res.RESC <- simresult$final.res.RESC
    POWsim1.RESC <- data.frame(final.res.RESC)


    POWsim1.RESC <- POWsim1.RESC %>%
      mutate(Effect=EFFECTRatio)


    POWsim1.RESC.long <- POWsim1.RESC %>%
      pivot_longer(-Effect, names_to = "Test", values_to = "Power")

    POWsim1.RESC.long$Test <- factor(POWsim1.RESC.long$Test, levels = unique(POWsim1.RESC.long$Test))


    ###

    POWsim1.ORIG.long$Score <- scoretype[1]
    POWsim1.RESC.long$Score <- scoretype[2]

    POWsim1.BOTH.long <- rbind(POWsim1.ORIG.long,POWsim1.RESC.long)

    POWsim1.BOTH.long$Test <- factor(POWsim1.BOTH.long$Test, levels = unique(POWsim1.BOTH.long$Test))


    ###


    # Define a function to create legend labels with subscripted values and EFFECTRatio values
    create_legend_labels <- function(values, ratios) {
      lapply(seq_along(values), function(i) {
        formatted_ratio <- sprintf("%.2f", ratios[i])
        bquote(rho[.(i)] == .(formatted_ratio))

      })
    }


    legend_labels <- create_legend_labels(1:length(EFFECTRatio), EFFECTRatio)

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



  }

}



