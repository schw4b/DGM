# Multiregression Dynamic Models (MDM)

# User Guide

# Developer Gruide

## Working with Github

### Check out master
    git clone https://github.com/schw4b/mdm.git
  
### Check out a release
    git checkout tags/v1.0
  
### Change to master
    git checkout origin/master
    
### Switch branches
    git checkout master
    git checkout develop

### Configure ssh key authentication
    git remote set-url origin git@github.com:schw4b/mdm.git

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
    
    
