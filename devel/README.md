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
    Unit: microseconds
                            expr      min        lq       mean    median        uq       max neval
    dlm.filt.rh(Yt, t(Ft), 0.93) 14161.16 14418.930 15197.0242 14719.324 14990.811 18437.171   100
     dlmFiltCpp(Yt, t(Ft), 0.93)   285.78   306.124   362.2461   327.562   421.552   457.031   100

40-fold speed improvements of the C++ function compared to R.

    Unit: milliseconds
                                   expr        min         lq       mean     median         uq       max neval
    exhaustive.search(myts, 3, cpp = F) 11356.4994 11462.7909 11597.2890 11559.1692 11651.8035 13080.311   100
             exhaustive.search(myts, 3)   186.8651   190.8783   193.6154   191.9987   193.3332   297.144   100

60-fold speed improvement compared to the native R implementation (cpp=F).

## Run benchmarks

    setwd('~/workspace')
    library(devtools)
    load_all('mdm')
    library(microbenchmark)
    data("utestdata")
    microbenchmark(exhaustive.search(myts,3,cpp=F),exhaustive.search(myts,3))
