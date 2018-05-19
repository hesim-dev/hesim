# ifndef INTEGRATE_H
# define INTEGRATE_H
#include <RcppArmadillo.h>

class Trapezoid{
public:
  double operator()(double a, double b, double y_a, double y_b) const {
    return (b - a)/2 * (y_a + y_b);
  }
};

template <class T>
class TrapezoidFun{
public:
  double operator()(T f, double a, double b) const {
    return (b - a)/2 * (f(a) + f(b));
  }
};

template <class T>
class SimpsonFun{
  public:
    double operator()(T f, double a, double b) const {
      return (b - a)/6 * (f(a) + 4 * f((a + b)/2) + f(b));
    }
};

template <typename Func, typename InputIt>
double composite(InputIt x_first, InputIt x_last, InputIt y_first, Func rule){
  double area = 0;
  InputIt y_it = y_first;
  for (InputIt x_it = x_first; x_it != (x_last - 1); ++x_it){
    area += rule(*x_it, *(x_it + 1), *y_it, *(y_it + 1));
    ++y_it;
  }
  return area;
}

template <typename Func1, typename Func2, typename InputIt>
double composite_fun(Func1 f, InputIt first, InputIt last, Func2 rule){
  double area = 0;
  for (InputIt it = first; it != (last - 1); ++it){
    area += rule(f, *it, *(it + 1));
  }
  return area;
}

template <typename Func1, typename Func2, typename InputIt>
std::vector<double> cumulative_composite_fun(Func1 f, InputIt first, InputIt last, Func2 rule){
  double area_i = 0;
  std::vector<double> area_vec(std::distance(first, last));
  area_vec[0] = 0;
  for (InputIt it = (first + 1); it != last; ++it){
    area_i += rule(f, *(it - 1), *it);
    area_vec[it - first] = area_i;
  }
  return area_vec;
}

template <typename Func, typename InputIt>
inline double trapzfun(Func f, InputIt first, InputIt last){
  TrapezoidFun<Func> trapz;
  return composite_fun(f, first, last, trapz);
}

template <typename InputIt>
inline double trapz(InputIt x_first, InputIt x_last, InputIt y_first){
  Trapezoid trapz;
  return composite(x_first, x_last, y_first, trapz);
}

template <typename Func, typename InputIt>
std::vector<double> cumtrapzfun(Func f, InputIt first, InputIt last){
  TrapezoidFun<Func> trapz;
  return cumulative_composite_fun(f, first, last, trapz);
}

template <typename Func, typename InputIt>
double simpsfun(Func f, InputIt first, InputIt last){
  SimpsonFun<Func> simps;
  return composite_fun(f, first, last, simps);
}

template <typename Func, typename InputIt>
std::vector<double> cumsimpsfun(Func f, InputIt first, InputIt last){
  SimpsonFun<Func> simps;
  return cumulative_composite_fun(f, first, last, simps);
}

# endif
