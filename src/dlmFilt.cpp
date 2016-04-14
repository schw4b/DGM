#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
  
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat dlmFiltCpp(NumericVector Yt_, NumericMatrix Ft_, int delta) {

  rowvec Yt(Yt_.begin(), Yt_.size(), false); // RcppArmadillo
  mat Ft(Ft_.begin(), Ft_.nrow(), Ft_.ncol(), false);
  
  int Nt = Yt.n_cols + 1; // the length of the time series + t0
  int  p = Ft.n_rows;     // the number of parents and one for an intercept (i.e. the number of thetas)
  
  rowvec m0 = zeros<rowvec>(p);
  mat CS0 = eye<mat>(p,p) * 3;
  float n0 = 0.001; float d0 = 0.001;
  
  rowvec Y(Yt);
  Y = join_rows(zeros<rowvec>(1), Yt);
  
  mat F1(Ft);
  F1 = join_horiz(zeros<colvec>(p), Ft);
  
  // Set up allocation matrices, including the priors
  mat mt = zeros<mat>(p,Nt);
  mt.col(0) = m0.t();
  
  cube Ct = zeros<cube>(p,p,Nt);
  Ct.slice(0) = CS0;
  cube CSt = zeros<cube>(p,p,Nt);
  
  cube Rt  = zeros<cube>(p,p,Nt);
  cube RSt = zeros<cube>(p,p,Nt);
  
  rowvec nt = zeros<rowvec>(Nt);
  nt(0) = n0;
  
  rowvec dt = zeros<rowvec>(Nt);
  dt(0) = d0;
  
  rowvec S = zeros<rowvec>(Nt);
  S(0) = dt[0]/nt[0];
  
  rowvec ft  = zeros<rowvec>(Nt);
  rowvec Qt  = zeros<rowvec>(Nt);
  rowvec ets = zeros<rowvec>(Nt);
  rowvec lpl = zeros<rowvec>(Nt);
  
  
  // Debug print to console
  Rcpp::Rcout << CS0;

  return S;
}
