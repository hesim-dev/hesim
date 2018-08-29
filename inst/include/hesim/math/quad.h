# ifndef HESIM_MATH_QUAD_H
# define HESIM_MATH_QUAD_H

# include <hesim/Rbase/integrate.h>

namespace hesim {

namespace math {

/** 
 * Internal details for hesim::math that should be ignored by external users.
 */
namespace detail {
/***************************************************************************//** 
 * Vectorize a functor.
 * Vectorize a functor so that it is computed over a @c array. Used in quad. 
 ******************************************************************************/ 		    
template <class Func>
class vectorize {
public:
  vectorize(Func f) :f_(f){}
  Func f_;
  void operator()(double *x, int n){
    for(int i = 0; i < n; i++){
      x[i] = f_(x[i]);
    }
  }
};	    

} // end namespace detail

/***************************************************************************//** 
 * Adaptive quadrature for a one-dimensional definite or indefinite integral.
 * Integrate a function using the same technique as in the @c R function
 * @c integrate from the @c stats package. Adapts the @c C functions from the 
 * @c stats package Rdqags and Rdqagi, which are based on the routines
 * @c dqags and @c dqagi from the @c Fortran library @c QUADPACK.  
 * @param[in] f A functor or lambda expression to integrate. The function must have
 * a single argument of type double.
 * @param[in] lower, upper The limits of integration. 
 * @param[out] abserr Estimate of the modulus of the absolute error.
 * @param[out] ier An integer equal to 0 if the routine terminated normally and
 * reliably; otherwise, the integer denotes a specific error message. See
 * the @c QUADPACK routines for more details.
 * @param[in] epsabs Absolute accuracy requested.
 * @param[in] epsrel Relative accuracy requested. 
 * @return The integral of the hazard function.
 ******************************************************************************/ 		    		    
template <typename Func>
double quad(Func f, double lower, double upper, double &abserr,
                 int &ier, double epsabs = 1e-6, double epsrel = 1e-6){
  double result;
  int neval;
  int limit = 100;
  int lenw = 4 * limit;
  int last;
  int *iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  double *work = (double *) R_alloc((size_t) lenw, sizeof(double));
  detail::vectorize<Func> vecf(f);
  if (std::isinf(lower) || std::isinf(upper)){ // start indefinite integral
    double bound = 0.0;
    int inf;
    if (!std::isinf(lower)){
      bound = lower;
      inf = 1;
    }
    else if (!std::isinf(upper)){
      bound = upper;
      inf = -1;
    }
    else{
      inf = 2;
    }
    Rbase::Rdqagi(vecf, &bound, &inf, &epsabs, &epsrel, &result,
                  &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  } // end indefinite integral
  else{ // start definite integral
    Rbase::Rdqags(vecf, &lower, &upper, &epsabs, &epsrel, &result,
                  &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  } // end definite integral
  return result;
}

} // end namespace math

} // end namespace hesim

# endif
