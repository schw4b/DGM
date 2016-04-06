context("Running MDM functions")

# ts 200x5 sub1 sim22
# utestdata benchmark values to test against

test_that("Running filtering distribution: 1 parent", {
  
  # # initial preparation of true values to test against later on
  # setwd("~/Data/NetSim/sim/smith-nodes/")
  # library(R.matlab)
  # d=readMat('sim22.mat')
  # myts=d$ts[1:200,]
  # rm(d)
  
  data("utestdata")
  Yt = myts[,1]
  expect_that(length(Yt), equals(200))
  
  Ft=array(1,dim=c(200,2))
  Ft[,2:ncol(Ft)]=myts[,2]
  
  # # initial preparation of true values to test against later on
  # a=dlm.filt(Yt,t(Ft),0.93)
  # utestdata=list()
  # utestdata$Np1.lpl=a$lpl
  # utestdata$Np1.lplsum=sum(a$lpl[15:200])
  # setwd("~/workspace/mdm/data/")
  # save(myts, utestdata, file = "utestdata.RData")

  a=dlm.filt(Yt,t(Ft),0.93)
  expect_that(length(a$lpl), equals(200))
  expect_that(sum(a$lpl[15:200]), equals(utestdata$Np1.lplsum))
  expect_that(a$lpl[15:200], equals(utestdata$Np1.lpl[15:200]))

})
