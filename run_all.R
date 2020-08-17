source("Ar_WCARC_pmf.R")
source("KB.R")
source("stepwise_rdf.R")
system("gfortran -O3 forceLMF.f")
system("./a.out")

source("LJ_overlap_sphere.R")
source("LMF_overlap_sphere.R")
source("compare_LMF_widom_simulation.R")