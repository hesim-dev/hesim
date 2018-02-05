// [[Rcpp::interfaces(r, cpp)]]
#include "distributions.h"
using namespace Rcpp;

/***************
* Free functions
***************/
double integrate_hazard(const Distribution * dist, double t){
    HazardFunc fun(dist);
    const double lower = 0, upper = t;
    double err_est; int err_code;
    return Numer::integrate(fun, lower, upper, err_est, err_code);
}

double quantile_numeric_work(const Distribution * dist, double p){
    QuantileNumericFunc func(dist, p);
    double lower = -1;
    double upper = 1;
    while(func(lower) * func(upper) >= 0){
        double interval = upper - lower;
        lower = lower - 0.5 * interval;
        upper = upper + 0.5 * interval;
    }
    double f_lower = func(lower);
    double f_upper = func(upper);
    double tol = 0.0001;
    int maxiter = 1000;
    return zeroin(lower, upper, f_lower, f_upper, func,
                  &tol, &maxiter);
}

double quantile_numeric(const Distribution * dist, double p){
  if ( p < 0 || p > 1){
    return NAN;
  }
  else if (p == 0) return R_NegInf;
  else if (p == 1) return R_PosInf;
  else{
    return quantile_numeric_work(dist, p);
  }
}

/**************************
* Exponential distribution
**************************/
Exponential::Exponential(double rate){
    rate_ = rate;
}

double Exponential::pdf(double x) const{
    return rate_ * exp(-rate_ * x);
}

double Exponential::cdf(double x) const{
    return 1 - exp(-rate_ * x); // R::pexp(x_, 1/rate_, 1, 0)
}

double Exponential::quantile(double p) const{
    return R::qexp(p, 1/rate_, 1, 0);
}

double Exponential::hazard(double x) const{
    return rate_;
}

double Exponential::cumhazard(double x) const{
    return rate_ * x;
}

double Exponential::random() const{
    return R::rexp(1/rate_);
}

/*********************
* Weibull distribution
*********************/
Weibull::Weibull(double shape, double scale){
    shape_ = shape;
    scale_ = scale;
}

double Weibull::pdf(double x) const{
    return R::dweibull(x, shape_, scale_, 0);
}

double Weibull::cdf(double x) const{
    return R::pweibull(x, shape_, scale_, 1, 0);
}

double Weibull::quantile(double p) const{
    return R::qweibull(p, shape_, scale_, 1, 0);
}

double Weibull::hazard(double x) const{
    return shape_ * pow(x/scale_, shape_ - 1)/scale_;
}

double Weibull::cumhazard(double x) const{
    return pow(x/scale_, shape_);
}

double Weibull::random() const{
    return R::rweibull(shape_, scale_);
}

/*******************
* Gamma distribution
*******************/
Gamma::Gamma(double shape, double rate){
    shape_ = shape;
    rate_ = rate;
}

double Gamma::pdf(double x) const{
    return R::dgamma(x, shape_, 1/rate_, 0);
}

double Gamma::cdf(double x) const{
    return R::pgamma(x, shape_, 1/rate_, 1, 0);
}

double Gamma::quantile(double p) const{
    return R::qgamma(p, shape_, 1/rate_, 1, 0);
}

double Gamma::hazard(double x) const{
    return Gamma::pdf(x)/(1 - Gamma::cdf(x));
}

double Gamma::cumhazard(double x) const{
    return -R::pgamma(x, shape_, 1/rate_, 0, 1);
}

double Gamma::random() const{
    return R::rgamma(shape_, 1/rate_);
}

/************************
* Lognormal distribution
***********************/
Lognormal::Lognormal(double meanlog, double sdlog){
    meanlog_ = meanlog;
    sdlog_ = sdlog;
}

double Lognormal::pdf(double x) const{
    return R::dlnorm(x, meanlog_, sdlog_, 0);
}

double Lognormal::cdf(double x) const{
    return R::plnorm(x, meanlog_, sdlog_, 1, 0);
}

double Lognormal::quantile(double p) const{
    return R::qlnorm(p, meanlog_, sdlog_, 1, 0);
}

double Lognormal::hazard(double x) const{
    return Lognormal::pdf(x)/(1 - Lognormal::cdf(x));
}

double Lognormal::cumhazard(double x) const{
    return -R::plnorm(x, meanlog_, sdlog_, 0, 1);
}

double Lognormal::random() const{
    return R::rlnorm(meanlog_, sdlog_);
}

/**********************
* Gompertz distribution
**********************/
// [[Rcpp::export]]
double qgompertz(double p, double shape, double rate) {
    double asymp = 1 - exp(rate/shape);
    if (shape == 0){
        return R::qexp(p, 1/rate, 1, 0);
    }
    else if (shape < 0 && p > asymp){
        return INFINITY;
    }
    else {
        return 1/shape * log(1 - shape * log(1 - p)/rate);
    }
}

// [[Rcpp::export]]
double rgompertz(double shape, double rate){
    double u = R::runif(0,1);
    return qgompertz(u, shape, rate);
}

Gompertz::Gompertz(double shape, double rate){
    shape_ = shape;
    rate_ = rate;
}

double Gompertz::pdf(double x) const{
    if (shape_ == 0){
        return R::dexp(x, 1/rate_, 0);
    }
    else{
        return rate_ * exp(shape_ * x) * exp(-rate_/shape_ * (exp(shape_ * x) -1));
    }
}

double Gompertz::cdf(double x) const{
    if (shape_ == 0){
        return R::pexp(x, 1/rate_, 1, 0);
    }
    else if (std::isinf(x)){
        return 1;
    }
    else{
        return 1 - exp(-rate_/shape_ * (exp(shape_ * x) - 1));
    }
}

double Gompertz::quantile(double p) const{
    return qgompertz(p, shape_, rate_);
}

double Gompertz::hazard(double x) const{
    return rate_ * exp(shape_ * x);
}

double Gompertz::cumhazard(double x) const{
    if (shape_ == 0){
        return rate_ * x;
    }
    else{
        return rate_/shape_ * expm1(shape_ * x);
    }
}

double Gompertz::random() const{
    return rgompertz(shape_, rate_);
}

/**************************
* Log-logistic distribution
**************************/
// [[Rcpp::export]]
double qllogis(double p, double shape, double scale, int lt = 1, int lg = 0){
    return exp(R::qlogis(p, log(scale), 1/shape, lt, lg));
}

// [[Rcpp::export]]
double rllogis(double shape, double scale){
    double u = R::runif(0,1);
    return qllogis(u, shape, scale);
}

LogLogistic::LogLogistic(double shape, double scale){
    shape_ = shape;
    scale_ = scale;
}

double LogLogistic::pdf(double x) const{
    return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    pow((1 + pow((x/scale_), shape_)), 2);
}

double LogLogistic::cdf(double x) const{
    return 1 - 1/(1 + pow(x/scale_, shape_));
}

double LogLogistic::quantile(double p) const{
    return exp(R::qlogis(p, log(scale_), 1/shape_, 1, 0));
}

double LogLogistic::hazard(double x) const{
    return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    (1 + pow((x/scale_), shape_));
}

double LogLogistic::cumhazard(double x) const{
    return -log(1 - LogLogistic::cdf(x));
}

double LogLogistic::random() const{
    return rllogis(shape_, scale_);
}

/********************************
* Generalized gamma distribution
********************************/
// [[Rcpp::export]]
double rgengamma(double mu, double sigma, double Q){
    if (Q == 0.0){
        return R::rlnorm(mu, sigma);
    }
    else{
        double w = log(pow(Q, 2) * R::rgamma(1/pow(Q, 2), 1))/Q;
        return exp(mu + sigma * w);
    }
}

GeneralizedGamma::GeneralizedGamma(double mu, double sigma, double Q){
    mu_ = mu;
    sigma_ = sigma;
    Q_ = Q;
}

double GeneralizedGamma::pdf(double x) const{
    if (Q_ != 0){
        double y = log(x);
        double w = (y - mu_)/sigma_;
        double Q2inv = 1/(Q_ * Q_);
        double logp = -log(sigma_ * x) + log(std::abs(Q_)) + Q2inv * log(Q2inv) +
        Q2inv * (Q_ * w - exp(Q_ * w)) - R::lgammafn(Q2inv);
        return exp(logp);
    } //
    else{
        return R::dlnorm(x, mu_, sigma_, 0);
    }
}

double GeneralizedGamma::cdf(double x) const{
    double y = log(x);
    double w = (y - mu_)/sigma_;
    double Q2inv = 1/(Q_ * Q_);
    double expnu = exp(Q_ * w) * Q2inv;
    if (Q_ > 0){
        return R::pgamma(expnu, Q2inv, 1, 1, 0);
    }
    else if (Q_ == 0){
        return R::plnorm(x, mu_, sigma_, 1, 0);
    }
    else{
        return 1 - R::pgamma(expnu, Q2inv, 1, 1, 0);
    }
}

double GeneralizedGamma::quantile(double p) const{
    if (Q_ == 0){
        return R::qlnorm(p, mu_, 1/(sigma_ * sigma_), 1, 0);
    }
    else {
        double gamma_quantile = R::qgamma(p, 1/(Q_ * Q_), 1, 1, 0);
        return exp(mu_ + sigma_ * (log(Q_ * Q_ * gamma_quantile)/Q_));
    }
}

double GeneralizedGamma::hazard(double x) const{
    return GeneralizedGamma::pdf(x)/(1 - GeneralizedGamma::cdf(x));
}

double GeneralizedGamma::cumhazard(double x) const{
    return -log(1 - GeneralizedGamma::cdf(x));
}

double GeneralizedGamma::random() const{
    return rgengamma(mu_, sigma_, Q_);
}

// [[Rcpp::export(name="C_rgengamma_vec")]]
std::vector<double> rgengamma_vec(int n, std::vector<double> mu,
                                  std::vector<double> sigma,
                                  std::vector<double> Q){
    std::vector<double> sample(n);
    int mu_size = mu.size();
    int sigma_size = sigma.size();
    int Q_size = Q.size();
    if (mu_size != sigma_size || mu_size != Q_size){
        Rcpp::stop("Length of mu, sigma, and Q must be the same");
    }
    for (int i = 0; i < n; ++i){
        int index = i % mu_size;
        GeneralizedGamma gengamma(mu[index], sigma[index], Q[index]);
        sample[i] = gengamma.random();
    }
    return(sample);
}

/************************
* Royston/Parmar Splines
*************************/
SurvSplines::SurvSplines(std::vector<double> gamma,
                         std::vector<double> knots,
                         std::string scale, std::string timescale){
    if (gamma.size() != knots.size()){
      Rcpp::stop("Length of gamma should equal number of knots.");
    }
    gamma_ = gamma;
    knots_ = knots;
    scale_ = scale;
    timescale_ = timescale;
    n_knots_ = knots.size();
    knot_max_ = *(knots.end() - 1);
    knot_min_ = *(knots.begin());
}

double SurvSplines::timescale_fun(double x) const{
    if (timescale_ == "log"){
        return log(x);
    }
    else if (timescale_ == "identity"){
        return x;
    }
    else{
        Rcpp::stop("Selected timescale is not available.");
    }
}

double SurvSplines::timescale_dx_fun(double x) const{
    if (timescale_ == "log"){
        return 1/x;
    }
    else if (timescale_ == "identity"){
        return 1;
    }
    else{
        Rcpp::stop("Selected timescale is not available.");
    }
}

double SurvSplines::basis_cube(double x) const{
    if (x <= 0) {
        return 0;
    }
    else {
        return x * x * x;
    }
}

double SurvSplines::basis_cube_dx(double x) const{
    if (x <= 0) {
        return 0;
    }
    else {
        return 3 * x * x;
    }
}

double SurvSplines::linear_predict(double x) const{
    double x_scaled = timescale_fun(x);
    std::vector<double> basis(n_knots_);
    basis[0] = 1; basis[1] = x_scaled;
    for (int j = 1; j < n_knots_ - 1; ++j){
        double lambda_j = (knot_max_ - knots_[j]) / (knot_max_ - knot_min_);
        basis[j + 1] = basis_cube(x_scaled - knots_[j]) - lambda_j *  basis_cube(x_scaled - knot_min_) -
        (1 - lambda_j) * basis_cube(x_scaled - knot_max_);
    }
    return std::inner_product(gamma_.begin(), gamma_.end(), basis.begin(), 0.0);
}

double SurvSplines::linear_predict_dx(double x) const {
    double x_scaled = timescale_fun(x);
    std::vector<double> basis_dx(n_knots_);
    basis_dx[0] = 0; basis_dx[1] = 1;
    for (int j = 1; j < n_knots_ - 1; ++j){
        double lambda_j = (knot_max_ - knots_[j]) / (knot_max_ - knot_min_);
        basis_dx[j + 1] = basis_cube_dx(x_scaled - knots_[j]) - lambda_j *  basis_cube_dx(x_scaled - knot_min_) -
        (1 - lambda_j) * basis_cube_dx(x_scaled - knot_max_);
    }
    return std::inner_product(gamma_.begin(), gamma_.end(), basis_dx.begin(), 0.0);
}

double SurvSplines::survival(double x) const{
    if (x <= 0){
        return 1; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard" | scale_ == "log_cumhazard"){
        return exp(-cumhazard(x));
    }
    else if (scale_ == "log_odds"){
        return 1/(1 + exp(linear_predict(x)));
    }
    else if (scale_ == "normal"){
        return R::pnorm(-linear_predict(x), 0, 1, 1, 0);
    }
    else {
        Rcpp::stop("Selected scale is not available.");
    }
}

double SurvSplines::hazard(double x) const{
    if (x <= 0){
        return 0; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard"){
        return exp(linear_predict(x));
    }
    else if (scale_ == "log_cumhazard"){
        return timescale_dx_fun(x) * linear_predict_dx(x) * exp(linear_predict(x));
    }
    else if (scale_ == "log_odds"){
        return timescale_dx_fun(x) * linear_predict_dx(x) * R::plogis(linear_predict(x), 0, 1, 1, 0);
    }
    else if (scale_ == "normal"){
        double lp = linear_predict(x);
        return timescale_dx_fun(x) * linear_predict_dx(x) * R::dnorm(-lp, 0, 1, 0)/R::pnorm(-lp, 0, 1, 1, 0);
    }
    else{
        Rcpp::stop("Selected scale is not available.");
    }
}

double SurvSplines::cumhazard(double x) const{
    if (x <= 0){
        return 0; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard"){
        return integrate_hazard(this, x);
    }
    else if (scale_ == "log_cumhazard"){
        return exp(linear_predict(x));
    }
    else if (scale_ == "log_odds"){
        return log1p(exp(linear_predict(x)));
    }
    else if (scale_ == "normal"){
        return -R::pnorm(-linear_predict(x), 0, 1, 1, 1);
    }
    else{
        Rcpp::stop("Selected scale is not available.");
    }
}

double SurvSplines::pdf(double x) const {
    if (x <= 0){
        return 0; // spline model is for time >= 0
    }
    double lp = linear_predict(x);
    double prob;
    if (scale_ == "log_hazard"){
        prob = survival(x) * hazard(x);
    }
    else if (scale_ == "log_cumhazard"){
        prob = timescale_dx_fun(x) * linear_predict_dx(x) * exp(lp - exp(lp));
    }
    else if (scale_ == "log_odds"){
        prob = timescale_dx_fun(x) * linear_predict_dx(x) * exp(lp - 2 * log(1 + exp(lp)));
    }
    else if (scale_ == "normal"){
        prob = timescale_dx_fun(x) * linear_predict_dx(x) * R::dnorm(lp, 0, 1, 0);
    }
    else{
        Rcpp::stop("Selected scale is not available.");
    }
    if (prob <= 0){
      prob = 0;
    }
    return prob;
}

double SurvSplines::cdf(double x) const{
    return 1 - survival(x);
}

double SurvSplines::quantile(double p) const{
  return quantile_numeric(this, p);
}

double SurvSplines::random() const{
    return quantile(R::runif(0, 1));
}

/************************
* Fractional polynomials
************************/
FracPoly::FracPoly(std::vector<double> gamma, std::vector<double> powers){
  gamma_ = gamma;
  powers_ = powers;
}

double FracPoly::basis_power(double x, double power) const{
  if (power == 0){
    return log(x);
  }
  else{
    return pow(x, power);
  }
}

std::vector<double> FracPoly::basis(double x) const{
  int n_powers = powers_.size();
  std::vector<double> basis(n_powers);
  basis[0] = basis_power(x, powers_[0]);
  double xp_old = basis[0];
  double xp_new;
  if (n_powers > 1){
    for (int i = 1; i < n_powers; ++i){
      if (powers_[i] == powers_[i - 1]){
        xp_new = log(x) * xp_old;
      }
      else {
        xp_new = basis_power(x, powers_[i]);
      }
      basis[i] = xp_new;
      xp_old = xp_new;
    }
  }
  return basis;
}

double FracPoly::linear_predict(double x) const{
  std::vector<double> b = basis(x);
  return std::inner_product(gamma_.begin(), gamma_.end(), b.begin(), 0.0);
}

double FracPoly::hazard(double x) const{
  return exp(linear_predict(x));
}

double FracPoly::cumhazard(double x) const{
  return integrate_hazard(this, x);
}

double FracPoly::cdf(double x) const{
  return 1 - exp(-cumhazard(x));
}

double FracPoly::pdf(double x) const{
  return hazard(x) * (1 - cdf(x));
}

double FracPoly::quantile(double p) const{
  return quantile_numeric(this, p);
}

double FracPoly::random() const{
  return quantile(R::runif(0, 1));
}
  
/*******************************
* Truncated normal distribution
********************************/
// [[Rcpp::export]]
double rtruncnorm(double mean, double sd, double lower, double upper){
    double  sample;
    sample = R::rnorm(mean, sd);
    while(sample < lower || sample > upper){
        sample = R::rnorm(mean, sd);
    }
    return sample;
}

/************************************
* Piecewise exponential distribution
*************************************/
// NOTE: rate in R::rexp is 1/rate in rexp!!!!!!!!
// [[Rcpp::export]]
double rpwexp (arma::rowvec rate, arma::rowvec time) {
    int T = rate.n_elem;
    double surv = 0.0;
    for (int t = 0; t < T; ++t){
        double rexp_t = R::rexp(1/rate(t));
        surv = time(t) + rexp_t;
        if (t < (T - 1)){
            if (surv < time(t + 1)){
                break;
            }
        }
    }
    return surv;
}

// Vectorized piecewise exponential
// [[Rcpp::export(name="C_rpwexp_vec")]]
std::vector<double> rpwexp_vec (int n, arma::mat rate, arma::rowvec time) {
    int b = rate.n_rows;
    std::vector<double> surv;
    surv.reserve(n);
    for (int i = 0; i < n; ++i){
        surv.push_back(rpwexp(rate.row(i % b), time));
    }
    return surv;
}

/*************************
* Categorical distribution
**************************/
// [[Rcpp::export]]
int rcat(arma::rowvec probs) {
    int k = probs.n_elem;
    double probs_sum = accu(probs);
    probs = probs/probs_sum;
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    int max = which_max(ans);
    return(max);
}

// [[Rcpp::export(name="C_rcat_vec")]]
arma::vec rcat_vec(int n, arma::mat probs){
    int b = probs.n_rows;
    arma::vec samp(n);
    for (int i = 0; i < n; ++i){
        samp(i) = rcat(probs.row(i % b));
    }
    return(samp);
}

/***********************
* Dirichlet distribution
************************/
// [[Rcpp::export]]
arma::rowvec rdirichlet(arma::rowvec alpha){
    int alpha_len = alpha.size();
    arma::rowvec x(alpha_len);
    for (int i = 0; i < alpha_len; ++i){
        x(i) = R::rgamma(alpha(i), 1);
    }
    return x/arma::sum(x);
}

// [[Rcpp::export(name="C_rdirichlet_mat")]]
arma::cube rdirichlet_mat(int n, arma::mat alpha){
    int J = alpha.n_rows;
    int K = alpha.n_cols;
    arma::cube samp(J, K, n);
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < J; ++j){
            samp.slice(i).row(j) = rdirichlet(alpha.row(j));
        }
    }
    return(samp);
}

/**************************************************
* Convert list of parameters from R to std::vector
**************************************************/
vecmats convert_distribution_parameters(std::string dist, Rcpp::List R_parlist){
    vecmats C_parlist;
    if (dist == "exponential" || dist == "exp"){
        C_parlist.push_back(as<arma::mat >(R_parlist["rate"]));
    }
    else if (dist == "weibull" || dist == "weibull.quiet" || dist == "llogis"){
        C_parlist.push_back(as<arma::mat >(R_parlist["shape"]));
        C_parlist.push_back(as<arma::mat >(R_parlist["scale"]));
    }
    else if (dist == "gompertz" || dist == "gamma"){
        C_parlist.push_back(as<arma::mat >(R_parlist["shape"]));
        C_parlist.push_back(as<arma::mat >(R_parlist["rate"]));
    }
    else if (dist == "lnorm"){
        C_parlist.push_back(as<arma::mat >(R_parlist["meanlog"]));
        C_parlist.push_back(as<arma::mat >(R_parlist["sdlog"]));
    }
    else if (dist == "gengamma"){
        C_parlist.push_back(as<arma::mat >(R_parlist["mu"]));
        C_parlist.push_back(as<arma::mat >(R_parlist["sigma"]));
        C_parlist.push_back(as<arma::mat >(R_parlist["Q"]));
    }
    else{
        Rcpp::stop("The selected distribution is not available.");
    }
    return C_parlist;
}

/*********************
* Select distribution
*********************/
Distribution * select_distribution(std::string dist_name,
                                   std::vector<double> parameters){
    Distribution *d;
    if (dist_name == "exponential" || dist_name == "exp"){
        d = new Exponential(exp(parameters[0]));
    }
    else if (dist_name == "weibull.quiet" || dist_name == "weibull"){
        d = new Weibull(exp(parameters[0]), exp(parameters[1]));
    }
    else if (dist_name == "gamma"){
        d = new Gamma(exp(parameters[0]), exp(parameters[1]));
    }
    else if (dist_name == "lnorm"){
        d = new Lognormal(parameters[0], exp(parameters[1]));
    }
    else if (dist_name == "gompertz"){
        d = new Gompertz(parameters[0], exp(parameters[1]));
    }
    else if (dist_name == "llogis"){
        d = new LogLogistic(exp(parameters[0]), exp(parameters[1])); 
    }
    else if (dist_name == "gengamma"){
        d = new GeneralizedGamma(parameters[0], exp(parameters[1]), parameters[2]); 
    }
    else{
        Rcpp::stop("The selected distribution is not available.");
    }
    return d;
}

/**************
* RCPP Modules
**************/
RCPP_MODULE(Distributions){
  class_<Distribution>("Distribution")
  .method("pdf", &Distribution::pdf)
  .method("cdf", &Distribution::cdf)
  .method("quantile", &Distribution::quantile)
  .method("hazard", &Distribution::hazard)
  .method("cumhazard", &Distribution::cumhazard)
  .method("random", &Distribution::random)
  ;

  class_<Exponential>("Exponential")
    .derives<Distribution>("Distribution")
    .constructor<double>()
    .method("pdf", &Exponential::pdf)
    .method("cdf", &Exponential::cdf)
    .method("quantile", &Exponential::quantile)
    .method("hazard", &Exponential::hazard)
    .method("cumhazard", &Exponential::cumhazard)
    .method("random", &Exponential::random)
  ;

  class_<Weibull>("Weibull")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Weibull::pdf)
    .method("cdf", &Weibull::cdf)
    .method("quantile", &Weibull::quantile)
    .method("hazard", &Weibull::hazard)
    .method("cumhazard", &Weibull::cumhazard)
    .method("random", &Weibull::random)
  ;
  
  class_<Gamma>("Gamma")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Gamma::pdf)
    .method("cdf", &Gamma::cdf)
    .method("quantile", &Gamma::quantile)
    .method("hazard", &Gamma::hazard)
    .method("cumhazard", &Gamma::cumhazard)
    .method("random", &Gamma::random)
  ;
  
  class_<Lognormal>("Lognormal")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Lognormal::pdf)
    .method("cdf", &Lognormal::cdf)
    .method("quantile", &Lognormal::quantile)
    .method("hazard", &Lognormal::hazard)
    .method("cumhazard", &Lognormal::cumhazard)
    .method("random", &Lognormal::random)
  ;
  
  class_<Gompertz>("Gompertz")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Gompertz::pdf)
    .method("cdf", &Gompertz::cdf)
    .method("quantile", &Gompertz::quantile)
    .method("hazard", &Gompertz::hazard)
    .method("cumhazard", &Gompertz::cumhazard)
    .method("random", &Gompertz::random)
  ;
  
  class_<LogLogistic>("LogLogistic")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &LogLogistic::pdf)
    .method("cdf", &LogLogistic::cdf)
    .method("quantile", &LogLogistic::quantile)
    .method("hazard", &LogLogistic::hazard)
    .method("cumhazard", &LogLogistic::cumhazard)
    .method("random", &LogLogistic::random)
  ;
  
  class_<GeneralizedGamma>("GeneralizedGamma")
    .derives<Distribution>("Distribution")
    .constructor<double, double, double>()
    .method("pdf", &GeneralizedGamma::pdf)
    .method("cdf", &GeneralizedGamma::cdf)
    .method("quantile", &GeneralizedGamma::quantile)
    .method("hazard", &GeneralizedGamma::hazard)
    .method("cumhazard", &GeneralizedGamma::cumhazard)
    .method("random", &GeneralizedGamma::random)
  ;
  
    class_<SurvSplines>("SurvSplines")
    .derives<Distribution>("Distribution")
    .constructor<std::vector<double>, std::vector<double>, std::string, std::string>()
    .method("linear_predict", &SurvSplines::linear_predict)
    .method("linear_predict_dx", &SurvSplines::linear_predict_dx)
    .method("pdf", &SurvSplines::pdf)
    .method("cdf", &SurvSplines::cdf)
    .method("quantile", &SurvSplines::quantile)
    .method("hazard", &SurvSplines::hazard)
    .method("cumhazard", &SurvSplines::cumhazard)
    .method("random", &SurvSplines::random)
  ;
    
     class_<FracPoly>("FracPoly")
    .derives<Distribution>("Distribution")
    .constructor<std::vector<double>, std::vector<double> >()
    .method("linear_predict", &FracPoly::linear_predict)
    .method("pdf", &FracPoly::pdf)
    .method("cdf", &FracPoly::cdf)
    .method("quantile", &FracPoly::quantile)
    .method("hazard", &FracPoly::hazard)
    .method("cumhazard", &FracPoly::cumhazard)
    .method("random", &FracPoly::random)
  ;
}
