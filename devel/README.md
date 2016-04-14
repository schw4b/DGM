# Multiregression Dynamic Models (MDM)
[![Build Status](https://travis-ci.org/schw4b/mdm.png?branch=develop)](https://travis-ci.org/schw4b/mdm)

# Developer Guide

## Working with Github

### Check out master
    git clone https://github.com/schw4b/mdm.git
  
### Check out a release
    git checkout tags/v1.0
  
### Change to master
    git checkout origin/master
    
### Switch branches / push to develop branch
    git checkout master
    git checkout develop
    git push origin develop

### Configure ssh key authentication
    git remote set-url origin git@github.com:schw4b/mdm.git
    
### Merge develop to master
    git checkout master
    git pull               # to update the state to the latest remote master state
    git merge develop      # to bring changes to local master from your develop branch
    git push origin master # push current HEAD to remote master branch

## Unit Testing
Test functions are written for the *testthat* package and can be found in the folder `tests`.
    
    # load MDM
    setwd('~/workspace')
    library(devtools)
    load_all('mdm')
    # Run unit tests
    library(testthat)
    test_dir('mdm/tests', reporter = 'Summary')
    
    testthat results ================================================================
    OK: 6 SKIPPED: 0 FAILED: 0
    
## Benchmarks
    > microbenchmark(dlm.filt.rh(Yt,t(Ft),0.93),dlm.filt(Yt,t(Ft),0.93),dlmFiltCpp(Yt,t(Ft),0.93))
    Unit: microseconds
                             expr       min         lq       mean    median        uq        max neval
     dlm.filt.rh(Yt, t(Ft), 0.93) 13833.870 13985.2575 14922.7911 14332.736 15036.854  21639.979   100
        dlm.filt(Yt, t(Ft), 0.93)  8121.926  8232.8760  9821.7963  8482.194  8667.776 106913.517   100
      dlmFiltCpp(Yt, t(Ft), 0.93)   286.036   330.7755   394.5441   410.480   426.150    632.925   100

