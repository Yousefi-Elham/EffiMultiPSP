## code to prepare `simPSPdata` dataset
library(readr)
simPSPdata <- read_csv("data-raw/simPSPdata.csv",col_types = list(col_integer()))

usethis::use_data(simPSPdata, overwrite = TRUE,compress = "xz")

