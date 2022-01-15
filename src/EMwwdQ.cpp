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
SEXP EMwwdQ(SEXP z, SEXP p, SEXP betal, SEXP tl, SEXP wwl, float sigma) {
  NumericMatrix zz(z); // N by g
  NumericMatrix pp(p); // N by g
  NumericVector bbetal(betal); // length g
  NumericMatrix ttl(tl); // N by Sl
  NumericVector wwwl(wwl); // length Sl for random effect l
  float sd(sigma);
  NumericVector temp(zz.nrow()); // length N
  NumericVector result(wwwl.size()); // length Sl for random effect l

  for(int i=0; i<zz.nrow(); i++){
    // i=1,...,N policyholders
    temp(i) = sum((zz(i,_)-pp(i,_))*bbetal); // sum over j's
  }
  for(int s=0; s<wwwl.size(); s++){
    // s=1,...,Sl clusters for l-th random effect
    result(s) = sum(ttl(_,s)*temp);
  }
  result = result-wwwl/(sd*sd); // divided by variance
  return result;
}
