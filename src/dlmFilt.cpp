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

  colvec Yt(Yt_.begin(), Yt_.size(), false); // RcppArmadillo
  mat Ft(Ft_.begin(), Ft_.nrow(), Ft_.ncol(), false);
  
  int Nt = Yt.n_rows + 1; // the length of the time series + t0
  int  p = Ft.n_rows;     // # the number of parents and one for an intercept (i.e. the number of thetas)
  
  //NumericVector Y=Yt;
  //Y.insert(Y.begin(),0);

  //mat F1(p,Nt);
  
  //Rcpp::Rcout << Nt;

  return Ft;
}
