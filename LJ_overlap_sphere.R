rm(list = ls())
library(tidyverse)
library(foreach)

sigA <- 0.34
sigW <- 0.316557
eps <- sqrt(0.99774*0.650194)
sig <- (sigA+sigW)/2
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

# r <- seq(0, 3, by = 0.01)
# u1r <- sapply(r, u1)
# plot(r, u1r)
d_KB <- read.table("d_KB.txt")$V1
a1 <- d_KB # the radius of sphere 1. this is the point where rdf_OA goes from 0 to non0
a2 <- 2.5*sig # range of interaction u2_AW

ds <- seq(0.3, 1.0, by = 0.01)

uLMF <- 
foreach(ids = 1:length(ds), .combine = "c") %do% {
  d <- ds[ids] # the distance between sphere 1 and 2
  
  r1 <- c(0, 0, 0) # the coordinate of sphere 1
  r2 <- c(0, 0, d) # the coordinate of sphere 2
  
  dr <- 0.005
  xpoints <- seq(-a1, a1, by = dr)
  
  
  npoints <- length(xpoints)
  allpoints <- matrix(0, nrow = npoints^3, ncol = 3)
  
  ic <- 1
  for (x in xpoints){
    for (y in xpoints){
      for (z in xpoints){
        allpoints[ic,] <- c(x,y,z)
        ic <- ic + 1
      }
    }
  }
  
  
  
  allpoints_to_r2 <- sweep(allpoints, 2, r2)
  
  dto1 <- sqrt(rowSums(allpoints**2))
  dto2 <- sqrt(rowSums(allpoints_to_r2**2))
  
  ind_overlap <- (dto1 < a1) & (dto2 < a2)
  
  points_overlaped <- allpoints[ind_overlap,]
  dto1_overlaped <- dto1[ind_overlap]
  dto2_overlaped <- dto2[ind_overlap]
  
  u1av <- mean(sapply(dto2_overlaped, u1))
  vol_overlaped <- nrow(points_overlaped)*dr**3
  
  rhobulk <- 7029/5.95**3 # water bulk density
  kT <- 2.5
  uwidom <- -2*u1av*vol_overlaped*rhobulk/kT
  uwidom
}

df_ruLMF <- data.frame(r = ds, u = uLMF)
saveRDS(df_ruLMF, file = "df_LMF_Ar_water")
