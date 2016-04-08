# Multiregression Dynamic Models (MDM)
[![Build Status](https://travis-ci.org/schw4b/mdm.png?branch=master)](https://travis-ci.org/schw4b/mdm)

# User Guide

# Developer Gruide

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
    
    
