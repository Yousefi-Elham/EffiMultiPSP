#' Read CSV duplicate
#' 
#' @param path path to file name
#' 
#' @return a \code{data.frame}
#' @export
#' @importFrom readr read_csv
#' 
#' @examples
#' CSV = system.file("extdata", "simPSPdata.csv", package= "effMultiend")
#' simPSPdata_read(CSV)

simPSPdata_read = function(path) {
  readr::read_csv(path)
}

