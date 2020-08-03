// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

// For accessing constant PI
#define _USE_MATH_DEFINES
#include <math.h>

class WeibullLogY: public Func
{
private:
  double k;
  double lambda;
public:
  WeibullLogY(double k_, double lambda_) : k(k_), lambda(lambda_) {}

  double operator()(const double& u) const
  {
    double temp = k/lambda * pow((u/lambda),k-1) * exp(-1.0*pow(u/lambda, k)) * log(u);
    // if(isinf(u)){
    //   temp = 0.0;
    // }else{
    //   temp = u * exp( log(k/lambda) - (k-1)*log(exp(u)/lambda) - pow(exp(u)/lambda, k) + u );
    // }
    if(isnan(temp) || isinf(temp)){
      return 0.0;
    }else{
      return temp;
    }
  }
};

Rcpp::List intWeibullLogY(double k, double lambda, double lower, double upper)
{
  WeibullLogY f(k, lambda);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code,
                               250, 1e-8, 1e-8, Integrator<double>::GaussKronrod71);
  return Rcpp::List::create(
    // Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

//' @export
// [[Rcpp::export]]
SEXP intWeibullLogYObs(SEXP k_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector k(k_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    result(j) = intWeibullLogY(k(0), lambda(0), lower(j), upper(j))[0];
  }

  return result;
}

//' @export
// [[Rcpp::export]]
SEXP intWeibullLogYLat(SEXP k_, SEXP lambda_, SEXP lower_, SEXP upper_)
{
  NumericVector k(k_);
  NumericVector lambda(lambda_);
  NumericVector lower(lower_);
  NumericVector upper(upper_);

  NumericVector result(lower.length());

  for(int j=0; j<lower.length(); j++){
    double temp1 = intWeibullLogY(k(0), lambda(0), R_NegInf, lower(j))[0];
    double temp2 = intWeibullLogY(k(0), lambda(0), upper(j), R_PosInf)[0];
    result(j) = temp1 + temp2;

  }

  return result;
}

