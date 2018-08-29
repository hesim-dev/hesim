# ifndef HESIM_CHECK_R_INFINITY_H
# define HESIM_CHECK_R_INFINITY_H

namespace hesim{

/**
 * @ingroup general
 * Check whether an object has the special value @c R_PosInf or @c R_NegInf. If
 * the value is @c R_PosInf, the value is changed to the @c C++ constant INFINITY and
 * if the value is @c R_NegInf, the value is changed to the @c C+ constant -INFINITY;
 * otherwise, the value is unchanged.
 * @param x An object to check
 * @return None
 */
inline void check_R_infinity(double &x){
  if (!R_FINITE(x)){
    if(x == R_PosInf){
      x = INFINITY;
    }
    else{
      x = -INFINITY;
    }
  }
}

} // end namespace hesim

# endif