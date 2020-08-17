rm(list = ls())
library(tidyverse)
library(foreach)

a1 <- 1 # the radius of sphere 1
a2 <- 1 # the radius of sphere 2
d <- 1.5 # the distance between sphere 1 and 2

r1 <- c(0, 0, 0) # the coordinate of sphere 1
r2 <- c(0, 0, d) # the coordinate of sphere 2

dr <- 0.01
xpoints <- seq(-a1, a1, by = dr)
# allpoints <- numeric()
# 
# for (x in xpoints){
#   for (y in xpoints){
#     for (z in xpoints){
#       allpoints <- rbind(allpoints, c(x,y,z))
#     }
#   }
# }

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

vol_overlaped <- nrow(points_overlaped)*dr**3

vol_analytic <- pi*(a1 + a2 - d)^2*(d^2 + 2*d*a2 - 3*a2*a2 + 2*d*a1 + 6*a1*a2 - 3*a1*a1)/(12*d)
