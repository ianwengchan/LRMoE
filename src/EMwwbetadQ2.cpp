#include <Rcpp.h>
using namespace Rcpp;

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
SEXP EMwwbetadQ2(SEXP zj, SEXP pj, SEXP p, double betajl, SEXP wbetal, SEXP wl, SEXP tl) {
  NumericVector zzj(zj); // length N
  NumericVector ppj(pj); // length N
  NumericMatrix pp(p);   // N by g matrix
  NumericVector weightedbeta(wbetal); // length N
  NumericVector wwl(wl); // length Sl
  NumericMatrix ttl(tl); // N by Sl matrix
  NumericVector result(ttl.ncol()); // length Sl

  for(int s=0; s<result.length(); s++){
    // s=1,...,Sl clusters of l-th random effect
    result(s) = sum((zzj - ppj - ppj*ttl(_,s)*wwl(s)*(betajl - weightedbeta))*ttl(_,s));
  }
  return result;
}
