// Calculate the Log Predictive Likelihood for a specified set of parents and a fixed delta.
// Copyright (C) 2016 Simon Schwab, Ruth Harbord, and Thomas Nichols.

// Usefull documentation:
// http://dirk.eddelbuettel.com/code/rcpp.armadillo.html

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::rowvec dlmFiltCpp(NumericVector Yt_, NumericMatrix Ft_, double delta, double m0_, double CS0_, double n0, double d0) {
  
  rowvec Yt(Yt_.begin(), Yt_.size(), false); // reuses memory and avoids extra copy
  mat Ft(Ft_.begin(), Ft_.nrow(), Ft_.ncol(), false);
  
  uword Nt = Yt.n_cols + 1; // the length of the time series + t0
  uword  p = Ft.n_rows;     // the number of parents and one for an intercept (i.e. the number of thetas)
  
  rowvec m0 = zeros<rowvec>(p);
  m0.ones();
  m0 = m0 * m0_;
  mat CS0 = eye<mat>(p,p) * CS0_;
  
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
  rowvec lpl = zeros<rowvec>(Nt-1);
  
  rowvec prod(1); // variables for result from matrix products
  rowvec QSt(1);
  
  double et;
  colvec At(p);
  
  // Updating
  
  for(uword t=1; t<Nt; ++t) { // unsigned signed problem
    
    // Posterior at {t-1}: (theta_{t-1}|D_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1} x d_{t-1}/n_{t-1}]
    // Prior at {t}: (theta_{t}|D_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, R_{t}]
    // D_{t-1} = y_{1},...,y_{t-1}
    
    // R*_{t} ~ C*_{t-1}/delta
    RSt.slice(t) = Ct.slice(t-1) / (S(t-1)*delta);
    Rt.slice(t) = RSt.slice(t) * S(t-1);
    // One-step forecast: (Y_{t}|D_{t-1}) ~ T_{n_{t-1}}[f_{t}, Q_{t}]
    prod = F1.col(t).t() * mt.col(t-1);
    ft(t) = prod(0);
    QSt = 1 + F1.col(t).t() * RSt.slice(t) * F1.col(t);
    Qt(t) = QSt(0) * S(t-1);
    et = Y(t) - ft(t);
    ets(t) = et / sqrt(Qt(t));
    
    // # Posterior at t: (theta_{t}|D_{t}) ~ T_{n_{t}}[m_{t}, C_{t}]
    // D_{t} = y_{1},...,y_{t}
    At = (RSt.slice(t) * F1.col(t))/QSt(0);
    mt.col(t) = mt.col(t-1) + (At*et);
    
    nt(t) = nt(t-1) + 1;
    dt(t) = dt(t-1) + (et*et)/QSt(0);
    S(t)=dt(t)/nt(t);
    
    CSt.slice(t) = RSt.slice(t) - (At * At.t())*QSt(0);
    Ct.slice(t) = S(t)*CSt.slice(t);
    
    // Log Predictive Likelihood 
    // in the original R code the asignment was to lpl(t) and the first values was discarded, so we use here lpl(t-1)
    lpl(t-1) = lgamma((nt(t-1)+1)/2)-lgamma(nt(t-1)/2)-0.5*log(PI*nt(t-1)*Qt(t))-((nt(t-1)+1)/2)*log(1+(1/nt(t-1))*(et*et)/Qt(t));
  }
  
  // Debug print to console
  // Rcout <<  Ct.slice(t) << '\n';
  
  return lpl;
}
