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
 * Maximum less than or equal to value.
 * Return the maximum element in a sorted vector less than or equal to @p value. Note that
 * the vector must be sorted from smallest to largest.
 * @param first, last Forward iterators defining the range to examine 
 * @param value The value that the maximum element must be less than.
 * @return Iterator pointing to the maximum element in the vector less than
 * or equal to @p value. 
 */
template <class ForwardIt, class T>
inline ForwardIt max_lt(ForwardIt first, ForwardIt last, const T& value){
  auto ub = std::upper_bound(first, last, value);
  return ub - 1;
}

/***************************************************************************//** 
 * Transition matrix.
 * A class for summarizing possible health state transitions in a multi-state 
 * model. 
 ******************************************************************************/
class trans_mat {
private:
  std::vector<std::vector<int> > trans_id_; ///< A vector of vectors. The outer vector 
  ///< denotes the starting state and each inner 
  ///< vector denotes the transition id
  ///< (indexed from 1 to patient::n_trans_)
  ///< corresponding to the possible transitions
  ///< from that state.
  
  std::vector<std::vector<int> > to_; ///< A vector of vectors. The outer vector
  ///< denotes the starting state and each inner
  ///< vector denotes a state that can be 
  ///< transitioned to.
  
  /** 
   * Count the number of non missing elements in the matrix. 
   * The number of non missing elements is equal to the number of 
   * possible transitions. 
   * @param m The same transition matrix as in the constructor.
   */                                          
  int count_non_nan(arma::mat m){
    int sum_non_nan = 0;
    for (int i = 0; i < m.n_rows; ++i){
      for (int j = 0; j < m.n_cols; ++j){
        if (!std::isnan(m(i, j))){
          ++sum_non_nan;
        }
      } // end loop over columns
    } // end loop over rows
    return sum_non_nan;
  }
  
  /** 
   * Determine whether each health state is absorbing.
   * @param trans A vector of vectors of the same format as trans_mat::trans_id_. Should
   * only be called after trans_ has been initialized. 
   */     
  std::vector<bool> is_absorbing(std::vector<std::vector<int> > trans){
    std::vector<bool> absorbing(trans.size());
    for (int i = 0; i < trans.size(); ++i){
      if(trans[i].size() > 0){
        absorbing[i] = false;
      }
      else {
        absorbing[i] = true;
      }
    } // end loop over states
    return absorbing;
  }
  
public:
  int n_trans_; ///< The total number of possible transitions.
  int n_states_; ///< The number of total health states.
  std::vector<bool> absorbing_; ///< A vector indicating whether each state is absorbing.
  ///< A state is absorbing if a row in the transition matrix
  ///< has all NAs. 
  
  /** 
   * The constructor.
   * @param m A matrix of integers indicating allowed transitions in a multi-state model
   *  in the format from the @c R package @c mstate. See
   * the argument "trans" in @c msprep in the @c mstate documentation.
   * @param R_index If TRUE, then transition ids in the matrix are assumed to be from R, 
   * and re-indexed to start from 0 (rather than 1).
   */                                       
  trans_mat(arma::mat m, bool R_index = true) {
    // Initialize n_trans_ and n_states_
    n_trans_ = count_non_nan(m);
    n_states_ = m.n_rows;
    
    // Initialize trans_ and to_
    for (int i = 0; i < m.n_rows; ++i){
      arma::rowvec m_row = m.row(i);
      std::vector<int> trans_i;
      std::vector<int> to_i;
      for (int j = 0; j < m_row.n_elem; ++j){
        if(!std::isnan(m_row(j))){
          if (R_index){
            trans_i.push_back(m_row(j) - 1); 
          }
          else{
            trans_i.push_back(m_row(j)); 
          }
          to_i.push_back(j);
        }
      } // end loop over columns
      to_.push_back(to_i);
      trans_id_.push_back(trans_i);
    } // end loop over rows
    
    // Initialize absorbing_
    absorbing_ = is_absorbing(trans_id_);
  }
  
  /** 
   * Return transition number ids. 
   * @param from_state The state to transition from.
   * @reurn A vector of the transitions numbers from the specified health state.
   */   
  std::vector<int> trans_id(int from_state) {
    return trans_id_[from_state];
  }
  
  /** 
   * Return states that can be transitioned to. 
   * @param from_state The state to transition from.
   * @return A vector of the transition states that can be transitioned to 
   * from the specified health state.
   */   
  std::vector<int> to(int from_state) {
    return to_[from_state];
  }  
  
};

  /***************************************************************************/
  /**
   * Integer bound for the x being in the vector range such that x<=range[1] -> 0
   * and range[i] < x <= range[i+1] -> i
   */
  inline int hesim_bound(double x, std::vector<double> range) {
    auto lower = std::lower_bound(range.begin(), range.end(), x);
    if (lower==range.end() || (*lower >= x && x > range[0]))
      lower--;
    return std::distance(range.begin(), lower);
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
