#include "integrate.h"

template <class T>
class Trapezoid{
  public:
    double operator()(T f, double a, double b) const {
      return (b - a)/2 * (f(a) + f(b));
    }
};

template <class T>
class Simpson{
  public:
    double operator()(T f, double a, double b) const {
      return (b - a)/6 * (f(a) + 4 * f((a + b)/2) + f(b));
    }
};

template <typename Func1, typename Func2>
double composite(Func1 f, std::vector<double> t, Func2 rule){
  int n = t.size();
  double area = 0;
  for (int i = 0; i < n - 1; ++i){
   area += rule(f, t[i], t[i + 1]);
 }
  return area;
}

template <typename Func1, typename Func2>
std::vector<double> cumulative_composite(Func1 f, std::vector<double> t, Func2 rule){
  int n = t.size(); 
  double area_i = 0;
  std::vector<double> area_vec(n);
  area_vec[0] = 0;
  for (int i = 1; i < n ; ++i){
   area_i += rule(f, t[i - 1], t[i]);
   area_vec[i] = area_i;
 }
  return area_vec;
}

template <typename Func>
double trapz(Func f, std::vector<double> t){
  Trapezoid<Func> trapz;
  return composite(f, t, trapz);
}

template <typename Func>
std::vector<double> cumtrapz(Func f, std::vector<double> t){
  Trapezoid<Func> trapz;
  return cumulative_composite(f, t, trapz);
}

template <typename Func>
double simps(Func f, std::vector<double> t){
  Simpson<Func> simps;
  return composite(f, t, simps);
}

template <typename Func>
std::vector<double> cumsimps(Func f, std::vector<double> t){
  Simpson<Func> simps;
  return cumulative_composite(f, t, simps);
}

class IntegrateTestFunc{
  public:
    double operator()(double x) const {
      return x * x;
    }
};

// [[Rcpp::export]]
double C_test_trapz(std::vector<double> t){
  IntegrateTestFunc fun;
  return trapz(fun, t);
}

// [[Rcpp::export]]
std::vector<double> C_test_cumtrapz(std::vector<double> t){
  IntegrateTestFunc fun;
  return cumtrapz(fun, t);
}

// [[Rcpp::export]]
double C_test_simps(std::vector<double> t){
  IntegrateTestFunc fun;
  return simps(fun, t);
}

// [[Rcpp::export]]
std::vector<double> C_test_cumsimps(std::vector<double> t){
  IntegrateTestFunc fun;
  return cumsimps(fun, t);
}