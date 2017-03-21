# multdyn: A package for Multiregression Dynamic Models (MDM)
[![Build Status](https://travis-ci.org/schw4b/multdyn.png?branch=master)](https://travis-ci.org/schw4b/multdyn)

Multiregression Dynamic Models (MDM) belong to the family of Dynamic Bayesian Networks. This package is designed to study  effective connectivity in functional MRI. MDM searches through all possible parent nodes for a specific node and provides an interpretable fit in terms of regression model for each network node (Costa et al., 2015).

Current research aims to fully characterize MDM using big data from the Human Connectome Project (HCP) in order to test validity and different aspects of reliability (test-retest reliability, out-of-sample reliability) in order to make this new method available for neuroimaging.

Reference to scientific articles see at the bottom.

## User Guide

### Install latest release from github
Run these commands only once, or if you want to upgrade to a new release.

    install.packages("devtools")
    library(devtools)
    install_github("schw4b/mdm@v1.3")

### Running a MDM example with simulated data
We load simulation data from Smith et al. (2011) of a 5-node network with 200 samples (time points) of one subject. Time series should be mean centered.

    library(mdmwarwick)
    data("utestdata")
    dim(myts)
    [1] 200   5

Now, let's do a full search across all possible parent models of 2^(n-1). Here, with n=5, we have 16 possible models.

    result=exhaustive.search(myts,3)
    result$model.store
              [,1]      [,2]      [,3]      [,4]      [,5]      [,6]     [,7]      [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]     [,15]     [,16]
    [1,]    1.0000    2.0000    3.0000    4.0000    5.0000    6.0000    7.000    8.0000    9.0000   10.0000   11.0000   12.0000   13.0000   14.0000   15.0000   16.0000
    [2,]    0.0000    1.0000    2.0000    4.0000    5.0000    1.0000    1.000    1.0000    2.0000    2.0000    4.0000    1.0000    1.0000    1.0000    2.0000    1.0000
    [3,]    0.0000    0.0000    0.0000    0.0000    0.0000    2.0000    4.000    5.0000    4.0000    5.0000    5.0000    2.0000    2.0000    4.0000    4.0000    2.0000
    [4,]    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.000    0.0000    0.0000    0.0000    0.0000    4.0000    5.0000    5.0000    5.0000    4.0000
    [5,]    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    5.0000
    [6,] -353.6734 -348.1476 -329.9119 -359.3472 -346.9735 -348.4044 -355.495 -347.4541 -337.3578 -333.6073 -353.8772 -349.6878 -346.9843 -358.2056 -341.7462 -355.5821
    [7,]    0.5000    0.7300    0.6800    0.7700    0.7600    0.7800    0.800    0.7900    0.7300    0.7700    0.8100    0.7900    0.8300    0.8400    0.7800    0.8300

The table colums are the 16 different models. First row indicates model number, rows 2-5 the parents, row 6 the log predictive likelihood (LPL), and row 7 the discount factor (delta). To get the winning model, we simply maximaze across LPL.

    which.max(result$model.store[6,])
    [1] 3

Model number 3 with node 2 as a parent is most likely.

### Analysis on the subject level
We do a full search on the subject level (exhautive search on each node). The list returned contains all the models, the winning models, and an adjacency matrix of the network.

    s=subject(myts)
    names(s)
    [1] "models" "winner" "adj"

### Plot network
    library(igraph)
    plot.net(s$adj)

![Network example](https://cloud.githubusercontent.com/assets/11832548/15656327/d1f13cf4-269d-11e6-95cb-c6bc1f190cf6.png)

## HPC guide (high performance computing)

### Timing estmates
Estimates for a time-series with 1200 samples (HCP), and for a 2.8GHz CPU.

| No. of nodes  | Time     |
| ------------- |:--------:|
| 3             | 0.2 sec  |
| 4             | 0.5 sec  |
| 5             | 1 sec  |
| 6             | 2.2 sec  |
| 7             | 4.8 sec  |
| 8             | 10 sec  |
| 9             | 22 sec  |
| 10            | 48 sec  |
| 11            | 1 min 46 sec  |
| 12            | 3 min 44 sec  |
| 13            | 8 min  8 sec  |
| 15            | 38 min |
| 20            | 30 hours |
| 25            | 58 days |

Timings are for one node only. To estimate the full network (all parents of all the nodes, the numbers above have to be multiplied by the number of nodes (e.g., a 10 node network takes approx. 10 minutes)


### Install or update  MDM on Buster super-computer (buster.stats.warwick.ac.uk)
From bash terminal on buster run:

    $ module load gcc
    $ module load R/3.2.4
    $ R

And from R run:

    install.packages("devtools") # run only once
    # type 22 and then 1 to select http mirror (https not working)
    library(devtools)
    install_github("schw4b/mdm@v1.4")
    library(mdmwarwick)

## References
1. Costa, L., Smith, J., Nichols, T., Cussens, J., Duff, E. P., and Makin, T. R. (2015). Searching Multiregression Dynamic Models of resting-state fMRI networks using integer programming. *Bayesian Analysis*, 10(2), 441–478. [doi:10.1214/14-BA913](http://dx.doi.org/10.1214/14-BA913).
2. Smith, S. M., Miller, K. L., Salimi-Khorshidi, G., Webster, M., Beckmann, C. F., Nichols, T. E., et al. (2011). Network modelling methods for FMRI. *NeuroImage*, 54(2), 875–91. [doi:10.1016/j.neuroimage.2010.08.063](http://dx.doi.org/10.1016/j.neuroimage.2010.08.063).
