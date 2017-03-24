# Multiregression Dynamic Models (MDM)
[![Build Status](https://travis-ci.org/schw4b/multdyn.png?branch=develop)](https://travis-ci.org/schw4b/multdyn)

# Developer Guide

## Working with Github

### Check out master
    git clone https://github.com/schw4b/multdyn.git

### Switch to develop branch
    git checkout develop

### Review, pull, commit and push
    git status
    git diff
    git pull
    git add foo
    git commit -m "Update"
    git push origin develop

### Configure ssh key authentication
    git remote set-url origin git@github.com:schw4b/multdyn.git

### Merge develop to master
    # update DESCRIPTION file with new version if applicable
    git checkout master
    git pull               # to update the state to the latest remote master state
    git merge develop      # to bring changes to local master from your develop branch
    git push origin master # push current HEAD to remote master branch

## Package building for CRAN
    make doc
    make build
    make file=multdyn_1.5.tar.gz check

## Unit Testing
Test functions are written for the *testthat* package and can be found in the folder `tests`.

    library(devtools)
    library(testthat)
    load_all('~/workspace/multdyn')
    # Run unit tests
    test_dir('~/workspace/multdyn/tests', reporter = 'Summary')

	testthat results ==============================================================================
	OK: 22 SKIPPED: 0 FAILED: 0

	DONE =========================================================================================

## Run benchmarks

    library(devtools)
    library(microbenchmark)
    load_all('~/workspace/multdyn')
    data("utestdata")
    microbenchmark(exhaustive.search(myts,3,cpp=F),exhaustive.search(myts,3))

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

## Magic stuff to resolve CRAN complaints

See https://github.com/RcppCore/Rcpp/issues/636

Download R-devel.tar.gz from https://cran.r-project.org/sources.html
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
