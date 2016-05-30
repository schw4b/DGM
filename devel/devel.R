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
install_github("schw4b/mdm@v1.1")

# @Source a C++ function ----
library(Rcpp)
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
microbenchmark(exhaustive.search(myts,3,cpp=F),exhaustive.search(myts,3), times=10L)

# 1200 volume bench
M=matrix(rnorm(1200*5), ncol=5)
microbenchmark(exhaustive.search(M,3,cpp=F),exhaustive.search(M,3), times=1L)

# @Quick and dirty benchmark ----
n = 1000;
start = Sys.time ()
#for (i in 1:n) {x = dlm.filt.rh(Yt,t(Ft),0.93)}  # 15.7 s
#for (i in 1:n) {x = dlm.filt(Yt,t(Ft),0.93)}     #  9.26 s
for (i in 1:n) {x = dlmFiltCpp(Yt,t(Ft),0.93)}    #  0.347 s
Sys.time () - start
# @Profiling ----
data("utestdata")
Rprof("profile.out")
#mymod = exhaustive.search(myts,2,cpp=F)
mymod = exhaustive.search(myts,2)
Rprof(NULL)
summaryRprof("profile.out")
