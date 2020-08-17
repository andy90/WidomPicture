rm(list = ls())
library(tidyverse)
df_pmf_diff <- readRDS("df_pmf_diff_simulation")

# the following are the results of the widom picture
deltaW_AW <- readRDS("df_LMF_Ar_water")
deltaW_WW <- readRDS("df_LMF_waterwater")

u_total <- deltaW_AW$u + deltaW_WW$u
deltaW <- data.frame(r = deltaW_AW$r, u = u_total)

# the following are the results of the LMF 
ruLMF_ArO <- read.table("ulmf_ArO.txt")
ruLMF_ArO_0 <- read.table("ulmf_ArO_0.txt")
ruLMF_OO <- read.table("ulmf_OO.txt")

kT <- 2.5
ruLMF_ArO$V2 <- ruLMF_ArO$V2/kT
ruLMF_ArO_0$V2 <- ruLMF_ArO_0$V2/kT
ruLMF_OO$V2 <- ruLMF_OO$V2/kT

r_all <- deltaW$r
uLMF_ArO <- approx(x = ruLMF_ArO$V1, y = ruLMF_ArO$V2, xout = r_all)$y
uLMF_ArO_0 <- approx(x = ruLMF_ArO_0$V1, y = ruLMF_ArO_0$V2, xout = r_all)$y
uLMF_OO <- approx(x = ruLMF_OO$V1, y = ruLMF_OO$V2, xout = r_all)$y

deltaW_AW_LMF <- data.frame(r = r_all, u = uLMF_ArO + uLMF_ArO_0)
deltaW_WW_LMF <- data.frame(r = r_all, u = uLMF_OO)
deltaW_LMF <- data.frame(r = r_all, u = uLMF_ArO + uLMF_ArO_0 + uLMF_OO)

deltaW_AW$type <- "solute-water vdw con. Widom"
deltaW_WW$type <- "water-water vdw con. Widom"
deltaW$type <- "total con. Widom"

deltaW_AW_LMF$type <- "solute-water vdw con. LMF"
deltaW_WW_LMF$type <- "water-water vdw con. LMF"
deltaW_LMF$type <- "total con. LMF"

df_pmf_diff$type <- "total con. simulation"

df_final <- do.call(rbind, list(deltaW_AW, deltaW_AW_LMF, deltaW_WW, deltaW_WW_LMF, deltaW, deltaW_LMF, df_pmf_diff))

df_final$type <- factor(df_final$type, levels = c("solute-water vdw con. Widom", "solute-water vdw con. LMF", "water-water vdw con. Widom", "water-water vdw con. LMF", "total con. Widom", "total con. LMF", "total con. simulation"))
p <- ggplot(data = df_final)
p <- p + geom_line(mapping = aes(x = r, y = u, color = type))
p <- p + geom_vline(xintercept = 0.35, linetype="dashed") 
p <- p + xlim(0.32, 1.0)
p <- p + ylim(-1.0, 1.1)
p <- p + xlab("r(nm)")
p <- p + ylab(TeX("$\\Delta \\omega(kT)"))
p

ggsave("widom_LMF_sim_ArAr.pdf")
