# Developer script to for unit testing, profiling, etc.
# Copyright (C) 2016 Simon Schwab, Ruth Harbord, and Thomas Nichols.

# @install packages ----
install.packages("microbenchmark")
install.packages("devtools")
install.packages("Rcpp")
install.packages("RcppArmadillo")
# @load package with devtools ----
setwd('~/workspace')
library(devtools)
load_all('mdm')
# unload("mdm/")
# @install directly from github ----
#install_github("schw4b/mdm")

# @Source a C++ function ----
sourceCpp("mdm/src/dlmFilt.cpp")
# @Run unit tests ----
library(testthat)
test_dir('mdm/tests', reporter = 'Summary')
# @C++ development ----
n=1;   # node to test
p=2:5; # parent node
nP=length(p);  # n parents

data("utestdata")
Yt = myts[,n]
Ft=array(1,dim=c(200,nP+1))
Ft[,2:ncol(Ft)]=myts[,p]

# R
x=dlm.filt(Yt,t(Ft),0.93)
# C++
sourceCpp('mdm/src/dlmFilt.cpp')
lpl=dlmFiltCpp(Yt,t(Ft),0.93)
boxplot(x$lpl - c(lpl[2:length(lpl)]))

# @Proper benchmark ----
library(microbenchmark)
# 4 parents
microbenchmark(dlm.filt.rh(Yt,t(Ft),0.93),dlmFiltCpp(Yt,t(Ft),0.93))

data("utestdata")
microbenchmark(exhaustive.search(myts,3,cpp=F),exhaustive.search(myts,3))

# @Quick and dirty benchmark ----
n = 1000;
start = Sys.time ()
#for (i in 1:n) {x = dlm.filt.rh(Yt,t(Ft),0.93)}  # 15.7 s
#for (i in 1:n) {x = dlm.filt(Yt,t(Ft),0.93)}     #  9.26 s
for (i in 1:n) {x = dlmFiltCpp(Yt,t(Ft),0.93)}    #  0.347 s
Sys.time () - start



microbenchmark(dlm.filt(Yt,t(Ft),0.93),       # 4 parents
               dlm.filt(Yt,t(Ft[,1:4]),0.93), # 3
               dlm.filt(Yt,t(Ft[,1:3]),0.93), # 2
               dlm.filt(Yt,t(Ft[,1:2]),0.93)) # 1 Effect of 1 to 5 nodes is 10%
# cmpfun     200  50 n of timepoints
# 4 parents  40%  41% reduction in time
# 1 parent   40%  44%
# @Profiling ----
data("utestdata")
Rprof("profile.out")
mymod = exhaustive.search(myts,2)
Rprof(NULL)
summaryRprof("profile.out")
