# 2016 Simon Schwab, s.schwab@warwick.ac.uk
# University of Warwick, United Kingdom
#
# developer script to for unit testing, profiling etc.

# install packages ----
install.packages("microbenchmark")
install.packages("devtools")
install.packages("Rcpp")
install.packages("RcppArmadillo")

# init ----
setwd('~/workspace')

library(devtools)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)

load_all('mdm')
# unload("mdm/")

# Call C++ function ----
sourceCpp("mdm/src/dlmFilt.cpp")
dlmFiltCpp(101)
# Run unit tests ----
library(testthat)
test_dir('mdm/tests', reporter = 'Summary')
# C++ development ----
n=1;   # node to test
p=2:5; # parent node
nP=length(p);  # n parents

data("utestdata")
Yt = myts[,n]
Ft=array(1,dim=c(200,nP+1))
Ft[,2:ncol(Ft)]=myts[,p]

x=dlm.filt.rh(Yt,t(Ft),0.93)

# C++
x=dlmFiltCpp(Yt,t(Ft),0.93)



# Benchmark ----
microbenchmark(dlm.filt.rh(Yt,t(Ft),0.93),dlm.filt(Yt,t(Ft),0.93))
microbenchmark(dlm.filt(Yt,t(Ft),0.93),       # 4 parents
               dlm.filt(Yt,t(Ft[,1:4]),0.93), # 3
               dlm.filt(Yt,t(Ft[,1:3]),0.93), # 2
               dlm.filt(Yt,t(Ft[,1:2]),0.93)) # 1 Effect of 1 to 5 nodes is 10%
# cmpfun     200  50 n of timepoints
# 4 parents  40%  41% reduction in time
# 1 parent   40%  44%
# Profiling ----
Rprof("profile.out")
mymod = exhaustive.search(ts,2)
Rprof(NULL)
summaryRprof("profile.out")
# Benchmarking speed ----
setwd("~/workspace/mdm/data/")
load('Sim22Sub1.rda')

mymod = exhaustive.search(ts,2)

# 11 seconds on a i5 1.80GHz
#  7.5 seconds -> 33% reduction
# crate data for unit testing ----
