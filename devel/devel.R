# 2016 Simon Schwab, s.schwab@warwick.ac.uk
# University of Warwick, United Kingdom

library(compiler)
source("~/workspace/mdm/R/mdm.R")

# Profiling ----
Rprof("profile.out")
mymod = exhaustive.search(ts,1)
Rprof(NULL)
summaryRprof("profile.out")

# Benchmarking speed ----
setwd("~/workspace/mdm/data/")
load('Sim22Sub1.rda')

mymod = exhaustive.search(ts,1)

# 11 seconds on a i5 1.80GHz
#  7.5 seconds -> 33% reduction