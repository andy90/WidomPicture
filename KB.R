# calculate the kirkwood buff integral

rm(list = ls())
library("latex2exp")
library("tidyverse")

rdfOA <- read.table("rdf_OA.xvg")

rdfOA_av <- 
  rdfOA %>%
  filter((V1 > 1.4) & (V1 < 2.0)) %>%
  summarise(bulk = mean(V2)) %>%
  pull()

rdfOA$V2 <- rdfOA$V2/rdfOA_av
KB <- 0
delr <- rdfOA$V1[2] - rdfOA$V1[1]
allr <- as.numeric()
allKB <- as.numeric()
for (i in 1:nrow(rdfOA)){
  r <- rdfOA$V1[i]+delr/2
  gr <- rdfOA$V2[i]
  KB <- KB + 4*pi*r**2*delr*(gr-1)
  allr <- c(allr, r)
  allKB <- c(allKB, KB)
}

df_KB <- data.frame(r = allr, KB = allKB)

KB_final <- 
df_KB %>%
  filter((r > 1.4) & (r < 2.0)) %>%
  summarise(bulk = mean(KB))
plot(allr, allKB)
d_KB <- (-KB_final/(4*pi/3))**(1.0/3.0)

write.table(d_KB, file = "d_KB.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
