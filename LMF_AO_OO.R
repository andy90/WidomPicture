rm(list = ls())
library(tidyverse)
library("latex2exp")


deltaW_AW <- readRDS("df_LMF_Ar_water")
deltaW_WW <- readRDS("df_LMF_waterwater")

u_total <- deltaW_AW$u + deltaW_WW$u
deltaW <- data.frame(r = deltaW_AW$r, u = u_total)

sigA <- 0.34
sigW <- 0.316557
eps <- 0.99774
sig <- sigA
u1 <- function(r){
  
  rwca <- sig*2.0**(1.0/6.0)
  
  u <- 0
  if (r <= rwca){
    u <- -eps
  }else{
    u <- 4*eps*((sig/r)**12-(sig/r)**6)
  }
  
  u
}

deltaw_sim <- -sapply(deltaW$r, u1)/2.5
df_sim <- data.frame(r = deltaW_AW$r, u = deltaw_sim, type = "sim")
deltaW$type <- "widom"
df_final <- rbind(deltaW, df_sim)

p <- ggplot(data = df_final)
p <- p + geom_line(mapping = aes(x = r, y = u, color = type))
p <- p + geom_vline(xintercept = 0.35, linetype="dashed") 
p <- p + xlim(0.32, 0.8)
p <- p + xlab("r(nm)")
p <- p + ylab(TeX("$\\Delta \\omega(kT)"))
p
