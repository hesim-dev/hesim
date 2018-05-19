#include <hesim/zeroin.h>

/*****************************************************
* Run an example from uniroot in the R stats library
*****************************************************/
class ZeroinTestFunc{
private:
  double a_;
public:
  ZeroinTestFunc(double a)
    :a_(a){}
  double operator()(double x) const {
    return x - a_;
  }
};

// [[Rcpp::export]]
double test_zeroin(){
  double lower = 0;
  double upper = 1;
  ZeroinTestFunc func(1.0/3.0);
  double f_lower = func(lower);
  double f_upper = func(upper);
  // double (*fcnPtr)(double) = zeroin_fun;
  double tol = 0.0001;
  int maxiter = 1000;
  double root = zeroin(lower, upper, f_lower, f_upper, func,
                   &tol, &maxiter);
  // return (double)maxiter;
  // return tol;
  return root;
}