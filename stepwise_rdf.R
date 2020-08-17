rm(list = ls())
library(tidyverse)

rdfOA <- read.table("rdf_OA.xvg")
d_KB <- read.table("d_KB.txt")$V1
rdfOA_stepwise <-
  rdfOA %>%
  mutate(g_new = if_else(V1 > d_KB, 1, 0)) %>%
  select(V1, g_new)

write.table(rdfOA_stepwise, "rdf_OA_stepwise.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
