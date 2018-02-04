# ifndef INTEGRATE_H
# define INTEGRATE_H
#include <RcppArmadillo.h>

template <typename Func>
double trapz(Func f, std::vector<double> t);

template <typename Func>
std::vector<double> cumtrapz(Func f, std::vector<double> t);

template <typename Func>
double simps(Func f, std::vector<double> t);

template <typename Func>
std::vector<double> cumsimps(Func f, std::vector<double> t);

# endif
