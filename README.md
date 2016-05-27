# Multiregression Dynamic Models (MDM)
[![Build Status](https://travis-ci.org/schw4b/mdm.png?branch=master)](https://travis-ci.org/schw4b/mdm)

Multiregression Dynamic Models (MDM) belong to the family of Dynamic Bayesian Networks and allows one to study effective connectivity in functional MRI. MDM searches through all possible parent nodes for a specific node and provides an interpretable fit in terms of regression model for each network node.

Current research aims to fully characterize MDM using big data from the Human Connectome Project (HCP) in order to test validity and different aspects of reliability (test-retest reliability, out-of-sample reliability).

## References
Costa, L., Smith, J., Nichols, T., Cussens, J., Duff, E. P., and Makin, T. R. (2015). Searching Multiregression Dynamic Models of resting-statetate fMRI networks using integer programming. Bayesian Analysis TBA, 1-38. [10.1214/14-BA913](http://dx.doi.org/10.1214/14-BA913).

## User Guide

### Install specific release from github

    > install.packages("devtools") # run only once
    > library(devtools)
    > install_github("schw4b/mdm@v1.1.1") # run only once
    > library(mdmwarwick)

### Find parents for node 3

    > data("utestdata")
    > result=exhaustive.search(myts,3)
    > result
    $model.store
              [,1]      [,2]      [,3]      [,4]      [,5]      [,6]     [,7]      [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]     [,15]     [,16]
    [1,]    1.0000    2.0000    3.0000    4.0000    5.0000    6.0000    7.000    8.0000    9.0000   10.0000   11.0000   12.0000   13.0000   14.0000   15.0000   16.0000
    [2,]    0.0000    1.0000    2.0000    4.0000    5.0000    1.0000    1.000    1.0000    2.0000    2.0000    4.0000    1.0000    1.0000    1.0000    2.0000    1.0000
    [3,]    0.0000    0.0000    0.0000    0.0000    0.0000    2.0000    4.000    5.0000    4.0000    5.0000    5.0000    2.0000    2.0000    4.0000    4.0000    2.0000
    [4,]    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.000    0.0000    0.0000    0.0000    0.0000    4.0000    5.0000    5.0000    5.0000    4.0000
    [5,]    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    5.0000
    [6,] -353.6734 -348.1476 -329.9119 -359.3472 -346.9735 -348.4044 -355.495 -347.4541 -337.3578 -333.6073 -353.8772 -349.6878 -346.9843 -358.2056 -341.7462 -355.5821
    [7,]    0.5000    0.7300    0.6800    0.7700    0.7600    0.7800    0.800    0.7900    0.7300    0.7700    0.8100    0.7900    0.8300    0.8400    0.7800    0.8300

    which.max(result$model.store[6,])

Model number 3 with node 2 as a parent is most likely.

  
### Remove package
    $ R --vanilla CMD REMOVE emov

### Install MDM on Buster Supercomputer
    # from terminal on buster
    $ module load gcc
    $ module load R/3.2.4
    $ R

    # within R
    > install.packages("devtools") # run only once
    # type 22 and then 1 to select http mirror (https not working)
    > library(devtools)
    > install_github("schw4b/mdm@v1.1.1") # run only once -- or if you want to update to a new release/tag
    > library(mdmwarwick)
