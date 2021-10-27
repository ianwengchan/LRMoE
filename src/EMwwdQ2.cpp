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
SEXP EMwwdQ2(SEXP p, SEXP betal, SEXP tl) {
  NumericMatrix pp(p); // N by g matrix
  NumericVector bbetal(betal); // length g
  NumericMatrix ttl(tl); // N by Sl matrix
  NumericVector prodsq(pp.nrow()); // length N, pp*beta
  NumericVector prod(pp.nrow()); // length N
  NumericMatrix temp(ttl.nrow(), ttl.ncol()); // N by Sl
  NumericMatrix result(ttl.ncol(), ttl.ncol()); // Sl by Sl matrix

  for(int i=0; i<pp.nrow(); i++){
    // i=1,...,N policyholders
    prodsq(i) = sum(pp(i,_)*bbetal); // sum over j's
    prod(i) = prodsq(i)*prodsq(i) - sum(pp(i,_)*bbetal*bbetal); // sum over j's
  }
  for(int j=0; j<result.nrow(); j++){
    // j=1,...,Sl clusters; not to mix up with j-th component
    temp(_,j) = ttl(_,j) * prod;
    for(int s=0; s<result.ncol(); s++){
      // s=1,...,Sl clusters
      result(j,s) = sum(temp(_,j)*ttl(_,s));
      if(s==j) result(j,s) = result(j,s) -1.0; // minus identity
    }
  }
  return result;
}
