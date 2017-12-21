# ifndef UTILS_H
# define UTILS_H
#include <RcppArmadillo.h>

typedef std::map<std::string, arma::mat> mapmat;
typedef std::map<std::string, std::vector<double> > mapvec;
typedef std::vector<arma::mat> vecmats;
vecmats list_to_vecmats(Rcpp::List l);
arma::mat matrix_byrow(arma::rowvec v, int nrow, int ncol);
arma::mat matrix_bycol(arma::rowvec v, int nrow, int ncol);

# endif
