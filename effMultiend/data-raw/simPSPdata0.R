## code to prepare `simPSPdata0` dataset
library(readr)
simPSPdata0 <- read_csv("data-raw/simPSPdata0.csv",col_types = list(col_integer()))

usethis::use_data(simPSPdata0, overwrite = TRUE,compress = "xz")

