rm(list = ls())
library("latex2exp")
library("tidyverse")
library("here")
pmf_2Ar <- read.table(here("ArgonWater/pmf_2Ar_water.xvg"))
pmf_2WCARC <- read.table(here("ArgonWater/pmf_2WCAAr_RC.xvg"))

pmf2Ar_av <- 
  pmf_2Ar %>%
  filter((V1 > 1.0) & (V1 < 1.1)) %>%
  summarise(bulk = mean(V2)) %>%
  pull()

pmf2WCARC_av <- 
  pmf_2WCARC %>%
  filter((V1 > 1.0) & (V1 < 1.1)) %>%
  summarise(bulk = mean(V2)) %>%
  pull()

r <- seq(0.3, 1.1, by = 0.005)
pmf_2Ar_std<- data.frame(approx(x = pmf_2Ar$V1, y = pmf_2Ar$V2, xout = r))
pmf_2WCARC_std<- data.frame(approx(x = pmf_2WCARC$V1, y = pmf_2WCARC$V2, xout = r))

pmf_2Ar_std$y <- pmf_2Ar_std$y - pmf2Ar_av
pmf_2WCARC_std$y <- pmf_2WCARC_std$y - pmf2WCARC_av

y_diff <- pmf_2Ar_std$y - pmf_2WCARC_std$y 
df_pmf_diff <- data.frame(r = r, u = y_diff)

saveRDS(df_pmf_diff, "df_pmf_diff_simulation")

