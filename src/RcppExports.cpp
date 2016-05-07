// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// markovCohortC
List markovCohortC(arma::cube P, arma::rowvec z0, int ncycles, arma::mat costs, arma::mat qol, double discount, arma::vec P_indx, arma::vec cost_indx, arma::vec qol_indx);
RcppExport SEXP cea_markovCohortC(SEXP PSEXP, SEXP z0SEXP, SEXP ncyclesSEXP, SEXP costsSEXP, SEXP qolSEXP, SEXP discountSEXP, SEXP P_indxSEXP, SEXP cost_indxSEXP, SEXP qol_indxSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::cube >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< int >::type ncycles(ncyclesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type costs(costsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type qol(qolSEXP);
    Rcpp::traits::input_parameter< double >::type discount(discountSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type P_indx(P_indxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cost_indx(cost_indxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type qol_indx(qol_indxSEXP);
    __result = Rcpp::wrap(markovCohortC(P, z0, ncycles, costs, qol, discount, P_indx, cost_indx, qol_indx));
    return __result;
END_RCPP
}
