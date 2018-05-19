# ifndef UTILS_H
# define UTILS_H
#include <RcppArmadillo.h>

typedef std::map<std::string, arma::mat> mapmat;
typedef std::map<std::string, std::vector<double> > mapvec;
typedef std::vector<double> vecdoubles;
typedef std::vector<vecdoubles> vecdoubles_2d;
typedef std::vector<std::string> vecstrings;
typedef std::vector<vecstrings> vecstrings_2d;
typedef std::vector<arma::mat> vecmats;
typedef std::vector<vecmats> vecmats_2d;
typedef std::vector<vecmats_2d> vecmats_3d;
typedef std::vector<arma::cube> vec_cubes;

vecmats list_to_vecmats(Rcpp::List l);
template <typename T> 
T list_to_vec(Rcpp::List l);
arma::mat matrix_byrow(arma::rowvec v, int nrow, int ncol);
arma::mat matrix_bycol(arma::rowvec v, int nrow, int ncol);

namespace Rcpp {
  template <> vecmats as(SEXP object);
  template <> vecmats_2d as(SEXP object);
  template <> vecmats_3d as(SEXP object);
  template <> vecstrings_2d as(SEXP object);
}

template <typename T>
void std_unique(std::vector<T> &vec){
  std::sort(vec.begin(), vec.end());
  vec.erase(unique(vec.begin(), vec.end()), vec.end());
}

template <typename T>
std::vector<T> add_constant(std::vector<T> &v, double value){
  transform(v.begin(), v.end(), v.begin(),
          bind2nd(std::plus<double>(), value)); 
  return v;
}

/****************
* Vector of zeros
****************/
template <typename T>
std::vector<double> zeros(int n){
  std::vector<T> v(n, 0.0);
  return v;
}

# endif
