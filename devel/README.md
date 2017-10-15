[![Build Status](https://travis-ci.org/schw4b/multdyn.png?branch=develop)](https://travis-ci.org/schw4b/multdyn)

# Developer Guide

## Manuals
- Rcpp: https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf
- RcppArmadillo: http://dirk.eddelbuettel.com/code/rcpp.armadillo.html

## Working with Github

### Check out master
    git clone https://github.com/schw4b/multdyn.git

### Configure ssh key authentication
    git remote set-url origin git@github.com:schw4b/multdyn.git

### Switch to develop branch
    git checkout develop

### Review, pull, commit and push
    git status
    git diff
    git pull
    git add foo
    git commit -m "Update"
    git push origin develop

### Merge develop to master
First, update DESCRIPTION file with new version

    git checkout master
    git pull               # to update the state to the latest remote master state
    git merge develop      # to bring changes to local master from your develop branch
    git push origin master # push current HEAD to remote master branch

## Package building for CRAN

### Required packages
    install.packages("devtools")
    install.packages("roxygen2")
    install.packages("data.table")
    install.packages("reshape2")
    install.packages("ggplot2")
    install.packages("RcppArmadillo")
    install.packages("testthat")


### Build a new package from the develop branch

Get `multdyn` and change to develop branch

    git clone https://github.com/schw4b/multdyn.git
    cd multdyn
    git checkout develop
    cd ..
    cp multdyn/devel/Makefile .

Build the package

    make clean
    cd multdyn; R -e 'devtools::document()'; cd ..
    make namespace
    make build
    make file=multdyn_1.6.1.tar.gz check

Install and test in R

    detach("package:multdyn", unload=TRUE)

    R CMD REMOVE multdyn
    R CMD INSTALL multdyn_1.6.1.tar.gz
    R

    library(multdyn)
    data(utestdata)
    result=exhaustive.search(myts,3)

## Unit Testing
Test functions are written for the *testthat* package and can be found in the folder `tests`. These tests do not have to be run manually as shown below, as they are run during the package building `make check`

    library(devtools)
    library(testthat)
    load_all('~/workspace/multdyn')
    # Run unit tests
    test_dir('~/workspace/multdyn/tests', reporter = 'Summary')

## Run benchmarks

    library(devtools)
    library(microbenchmark)
    load_all('~/workspace/multdyn')
    data("utestdata")
    microbenchmark(exhaustive.search(myts,3,cpp=F),exhaustive.search(myts,3), times=1)

Output:

     Unit: milliseconds
                                expr       min        lq      mean    median
    exhaustive.search(myts, 3, cpp = F) 7405.6748 7405.6748 7405.6748 7405.6748
             exhaustive.search(myts, 3)  238.3177  238.3177  238.3177  238.3177

30-fold speed improvement of the C++ function compared to R.

## Travis CI
The new version of Rcpp requires a higher version of g++. This can bes specified in `.travis.yml`. These lines will switch to a higher Ubuntu version with a newer g++/gcc version:

    # Ubuntu 14.04 Trusty support
    sudo: required
    dist: trusty

## Resolve CRAN complaints (obsolete)

As with the new `Rcpp` version this does not seem to be necessary anymore and I have removed it.

Following problem occured: https://github.com/RcppCore/Rcpp/issues/636

Download R-devel.tar.gz from https://cran.r-project.org/sources.html
```
mkdir ~/tmp
cd ~/tmp
wget https://stat.ethz.ch/R/daily/R-devel.tar.gz
tar -xvf R-devel.tar.gz
```
Then, in R, do the following:
```
source("~/tmp/R-devel/src/library/tools/R/sotools.R")
source("~/tmp/R-devel/src/library/tools/R/utils.R")
setwd("~/workspace/multdyn")

package_native_routine_registration_skeleton(".")
```
The function will output some code, copy and paste this to `src/register.c`
```
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP multdyn_dlmLplCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"multdyn_dlmLplCpp", (DL_FUNC) &multdyn_dlmLplCpp, 7},
    {NULL, NULL, 0}
};

void R_init_multdyn(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
```

Finally, this needs to go into `NAMESPACE`, see Makefile.
```
useDynLib(packagename, .registration = TRUE).
```
