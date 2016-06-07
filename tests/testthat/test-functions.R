context("Running MDM functions")

# ts 200x5 sub1 sim22
# utestdata benchmark values to test against

test_that("Running filtering distribution: 1 parent", {
  
  n=1;   # node to test
  p=2;   # parent node(s)
  nP=length(p);  # n parents
  
  # # initial preparation of true values to test against later on
  # setwd("~/Data/NetSim/sim/smith-nodes/")
  # library(R.matlab)
  # d=readMat('sim22.mat')
  # myts=d$ts[1:200,]
  # rm(d)
  
  data("utestdata")
  Yt = myts[,n]
  Ft=array(1,dim=c(200,nP+1))
  Ft[,2:ncol(Ft)]=myts[,p]
  
  # # initial preparation of true values to test against later on
  # a=dlm.filt(Yt,t(Ft),0.93)
  # utestdata=list()
  # utestdata$Np1.lpl=a$lpl
  # utestdata$Np1.lplsum=sum(a$lpl[15:200])
  # setwd("~/workspace/mdm/data/")
  # save(myts, utestdata, file = "utestdata.RData")
  
  a=dlm.filt.rh(Yt,t(Ft),0.93)
  expect_that(sum(a$lpl[15:200]), equals(utestdata$Np1.lplsum))
  expect_that(a$lpl, equals(utestdata$Np1.lpl))
  
  lpl=c(dlmFiltCpp(Yt,t(Ft),0.93))
  expect_equal(sum(lpl[15:200]), utestdata$Np1.lplsum)
  expect_equal(lpl, utestdata$Np1.lpl)
  
})

test_that("Running filtering distribution: 2 parents", {
  
  # sample(1:5,5, replace=F) # which node should i test?
  n=3;   # node to test
  p=c(1,4); # parent node
  nP=length(p);  # n parents
  
  data("utestdata")
  Yt = myts[,n]
  Ft=array(1,dim=c(200,nP+1))
  Ft[,2:ncol(Ft)]=myts[,p]
  
  # # initial preparation of true values to test against later on
  # utestdata$Np2.lpl=a$lpl
  # utestdata$Np2.lplsum=sum(a$lpl[15:200])
  # setwd("mdm/data/")
  # save(myts, utestdata, file = "utestdata.RData")
  
  a=dlm.filt.rh(Yt,t(Ft),0.93)
  expect_that(sum(a$lpl[15:200]), equals(utestdata$Np2.lplsum))
  expect_that(a$lpl, equals(utestdata$Np2.lpl))
  
  lpl=c(dlmFiltCpp(Yt,t(Ft),0.93))
  expect_equal(sum(lpl[15:200]), utestdata$Np2.lplsum)
  expect_equal(lpl, utestdata$Np2.lpl)
  
})

test_that("Running filtering distribution: 3 parents", {
  
  # sample(1:5,5, replace=F) # which node should i test?
  n=4;   # node to test
  p=c(1,2,3); # parent node
  nP=length(p);  # n parents
  
  data("utestdata")
  Yt = myts[,n]
  Ft=array(1,dim=c(200,nP+1))
  Ft[,2:ncol(Ft)]=myts[,p]
  
  # # initial preparation of true values to test against later on
  # utestdata$Np3.lpl=a$lpl
  # utestdata$Np3.lplsum=sum(a$lpl[15:200])
  # setwd("mdm/data/")
  # save(myts, utestdata, file = "utestdata.RData")
  
  a=dlm.filt.rh(Yt,t(Ft),0.93)
  expect_that(sum(a$lpl[15:200]), equals(utestdata$Np3.lplsum))
  expect_that(a$lpl, equals(utestdata$Np3.lpl))
  
  lpl=c(dlmFiltCpp(Yt,t(Ft),0.93))
  expect_equal(sum(lpl[15:200]), utestdata$Np3.lplsum)
  expect_equal(lpl, utestdata$Np3.lpl)
  
})

test_that("Running filtering distribution: 4 parents", {
  
  n=1;   # node to test
  p=2:5; # parent node
  nP=length(p);  # n parents
  
  data("utestdata")
  Yt = myts[,n]
  Ft=array(1,dim=c(200,nP+1))
  Ft[,2:ncol(Ft)]=myts[,p]
  
  # # initial preparation of true values to test against later on
  # utestdata$Np4.lpl=a$lpl
  # utestdata$Np4.lplsum=sum(a$lpl[15:200])
  # setwd("mdm/data/")
  # save(myts, utestdata, file = "utestdata.RData")
  
  a=dlm.filt.rh(Yt,t(Ft),0.93)
  expect_that(sum(a$lpl[15:200]), equals(utestdata$Np4.lplsum))
  expect_that(a$lpl, equals(utestdata$Np4.lpl))
  
  lpl=c(dlmFiltCpp(Yt,t(Ft),0.93))
  expect_equal(sum(lpl[15:200]), utestdata$Np4.lplsum)
  expect_equal(lpl, utestdata$Np4.lpl)
  
})

test_that("Exhaustive search, 5 node network", {
  
  data("utestdata")
  
  # # generate test values with original function for all 5 nodes
  # models = array(NA, dim=c(7,16,5))
  # for (n in 1:5) {
  #   mymod = exhaustive.search(myts,n)
  #   models[,,n] = mymod$model.store
  # }
  # utestdata$models=models
  # setwd("mdm/data/")
  # save(myts, utestdata, file = "utestdata.RData")
  
  # calculate every parent and compare
  for (n in 1:5) {
    mymod = exhaustive.search(myts,n)
    expect_equivalent(mymod$model.store,utestdata$models[,,n])
  }
})

test_that("center", {
  
  X = array(c(11,22,28,44,55,99), dim=c(2,2))
  M = X
  M[,1] = X[,1] - mean(X[,1])
  M[,2] = X[,2] - mean(X[,2])
  
  expect_equal(center(X), M)
})
