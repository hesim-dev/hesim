# ifndef HESIM_MATH_RIEMANN_H
# define HESIM_MATH_RIEMANN_H
#include <vector>

namespace hesim {

/** @ingroup math 
 * Functions and classes for numerical computing including integration 
 * and root finding.
 */
namespace math {

/***************************************************************************//** 
 * Cumulative Riemann sum.
 * Compute (midpoint) cumulative Riemann sum given function @c fun at the points in x.  
 * @param x_first, x_last Iterators defining the values of the points x.  
 * @param f Function used to compute sum.
 * @return The cumulative integral of the function approximated by the Riemann sum.
 ******************************************************************************/ 
template <typename InputIt, typename Func>
inline std::vector<double> cum_riemann(InputIt x_first, InputIt x_last, Func f){
  int n = std::distance(x_first, x_last);
  std::vector<double> sum(n);
  sum[0] = 0;
  double sum_i = 0;
  for (InputIt x_it = (x_first + 1); x_it != x_last; ++x_it){
    double step = *x_it - *(x_it - 1);
    double mid = *(x_it - 1) + step/2;
    sum_i += step * f(mid);
    sum[x_it - x_first] = sum_i; 
  }
  return sum;
}

} // end namespace math

} // end namespace hesim

# endif
