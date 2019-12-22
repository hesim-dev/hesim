# ifndef HESIM_UTILS_H
# define HESIM_UTILS_H
#include <RcppArmadillo.h>

namespace hesim{

typedef std::vector<arma::mat> vecmats;
typedef std::vector<vecmats> vecmats_2d;
typedef std::vector<vecmats_2d> vecmats_3d;
typedef std::vector<std::string> vecstrings;
typedef std::vector<vecstrings> vecstrings_2d;
typedef std::vector<arma::cube> vec_cubes;

/** 
 * @ingroup general
 * Internal details for hesim that should be ignored by external users.*/
namespace detail {

/**
 * Convert an Rcpp::List to a vector  
 * @param l An Rcpp::List 
 * @tparam T1 The class of the object to create (e.g., std::vector<arma::mat>).
 * @tparam T2 The class of each element of the vector (e.g., arma::mat).
 * @return A vector of type T1
 */
template <typename T1, typename T2> 
T1 list_to_vec(Rcpp::List l){
  T1 v;
  int n = l.size();
  v.reserve(n);
  for (int i = 0; i < n; ++i){
    v.push_back(Rcpp::as<T2 >(l[i]));
  }
  return v;
}

} //end namespace detail

/**
 * @ingroup general
 * Find the position of the largest element in the range
 * [first,last) from a container in the Standard 
 * Library.
 * @param first, last Forward iterators defining the range to examine 
 */
template <typename InputIt>
inline int max_element_pos(InputIt first, InputIt last) {
  auto it = std::max_element(first, last);
  return std::distance(first, it);
}

/**
 * @ingroup general
 * Sort a vector in the standard library and
 * erase duplicates.
 * @param v A vector.
 * @return None. 
 */
template <typename T>
inline void unique(std::vector<T> &v){
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

/**
 * @ingroup general
 * Add a constant value to a vector in the Standard
 * Library. 
 * @param v A vector. Should be of type integer
 * or double.
 * @param value A value to add to each element
 * in the vector. 
 * @return None.
 */
template <typename T>
inline void add_constant(std::vector<T> &v, double value){
  auto func = [value](T& d){ 
    return d+=value;
  };  
  std::for_each(v.begin(), v.end(), func);
}

/**
 * @ingroup general
 * Group counter for a sorted vector.
 * @param v A vector to group by. @p v MUST have been previously sorted.
 * @return An integer vector equal to 0 for the first group, 1 for the second group, 
 * and so on.
 */
template<typename T>
inline std::vector<int> grp_counter(std::vector<T> v){
  std::vector<int> v2(v.size());
  int grp = 0;
  v2[0] = grp;
  for (int i = 1; i < v.size(); ++i){
    if(v[i] != v[i-1]){
      ++grp;
    }
    v2[i] = grp; 
  }
  return v2;
}

/**
 * @ingroup general
 * Compute present value in continuous time. 
 * @param z A fixed value to be discounted.
 * @param r The discount rate
 * @param t1 Time at the start of the interval.
 * @param t2 Time at the end of the interval.
 * @return Present value with continuous compounding.  
 */
inline double pv(double z, double r, double t1, double t2){
  if (r == 0.0){
    return z * (t2 - t1);
  }
  else{
    return z * ((exp(-r * t1) - exp(-r * t2))/r); 
  }
}

/**
 * @ingroup general
 * Generate a sequence of numbers.
 * Generate a sequence of numbers starting at @p from and ending at @p to with a
 * step size of @p by. Note that caution should be taken because of possible 
 * floating point arithmetic errors.
 * @param from, to The starting and (maximal) end values of the sequence.
 * @param by Step size of the sequence. 
 * @return A sequence of numbers.
 */
inline std::vector<double> seq(double from, double to, double by){
  if ((from < to && by < 0) || (from > to && by > 0)){
    Rcpp::stop("Wrong sign in 'by' argument.");
  } 
  int size = int((to - from)/by) + 1;
  std::vector<double> result(size);
  result[0] = from;
  if (size > 1){
    for (int i = 1; i < size; ++i){
      result[i] = result[i - 1] + by;
    }
  }
  return result;
};

/**
 * @ingroup general
 * Check if named element in a list is NULL.
 * The named element is considered NULL if it is missing from the list
 * or it exists in the list but has a value of NULL. 
 * @param L The list.
 * @param name The name of the element to check.
 * @return true if the named element is NULL and false otherwise. 
 */
inline bool is_null(Rcpp::List L, const char * name){
  return  !L.containsElementNamed(name) || Rf_isNull(L[name]);
}

/**
 * @ingroup general
 * Maximum less than value.
 * Return the maximum element in a sorted vector less than @p value. Note that
 * the vector must be sorted from smallest to largest.
 * @param first, last Forward iterators defining the range to examine 
 * @param value The value that the maximum element must be less than.
 * @return Iterator pointing to the maximum element in the vector less than
 * @p value. 
 */
template <class ForwardIt, class T>
inline ForwardIt max_lt(ForwardIt first, ForwardIt last, const T& value){
  auto lb = lower_bound(first, last, value);
  if (lb == first){
    Rcpp::stop("There is no element in the vector less than 'value'");
  } 
  else{
    return lb - 1;
  }
}

} // end hesim namespace


/****************************
* Custom Rcpp::as converters
****************************/
namespace Rcpp {
  inline hesim::vecmats as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::detail::list_to_vec<hesim::vecmats, arma::mat>(l);
  }
  
  template <> inline hesim::vecmats_2d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::detail::list_to_vec<hesim::vecmats_2d, hesim::vecmats>(l);
  }
  
  template <> inline hesim::vecmats_3d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::detail::list_to_vec<hesim::vecmats_3d, hesim::vecmats_2d> (l);
  }
  
  template <> inline hesim::vecstrings_2d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::detail::list_to_vec<hesim::vecstrings_2d, hesim::vecstrings> (l);
  }
} // end Rcpp namespace



# endif
