// [[Rcpp::interfaces(r, cpp)]]
#include <hesim/distributions.h>

/**************************
* Exponential distribution
**************************/
hesim::Exponential::Exponential(double rate){
    rate_ = rate;
}

void hesim::Exponential::set_params(std::vector<double> params){
  rate_ = exp(params[0]);
}

double hesim::Exponential::pdf(double x) const{
    return rate_ * exp(-rate_ * x);
}

double hesim::Exponential::cdf(double x) const{
    return 1 - exp(-rate_ * x); // R::pexp(x_, 1/rate_, 1, 0)
}

double hesim::Exponential::quantile(double p) const{
    return R::qexp(p, 1/rate_, 1, 0);
}

double hesim::Exponential::hazard(double x) const{
    return rate_;
}

double hesim::Exponential::cumhazard(double x) const{
    return rate_ * x;
}

double hesim::Exponential::random() const{
    return R::rexp(1/rate_);
}

/*********************
* Weibull distribution
*********************/
hesim::Weibull::Weibull(double shape, double scale){
    shape_ = shape;
    scale_ = scale;
}

hesim::Weibull hesim::Weibull::create_from_Nma(double a0, double a1){
  double shape = a1 + 1;
  double scalePH = exp(a0)/shape;
  double scale = pow(scalePH, -1/shape);
  return Weibull(shape, scale);
}

void hesim::Weibull::set_params(std::vector<double> params){
  shape_ = exp(params[0]);
  scale_ = exp(params[1]);
}

double hesim::Weibull::pdf(double x) const{
    return R::dweibull(x, shape_, scale_, 0);
}

double hesim::Weibull::cdf(double x) const{
    return R::pweibull(x, shape_, scale_, 1, 0);
}

double hesim::Weibull::quantile(double p) const{
    return R::qweibull(p, shape_, scale_, 1, 0);
}

double hesim::Weibull::hazard(double x) const{
    return shape_ * pow(x/scale_, shape_ - 1)/scale_;
}

double hesim::Weibull::cumhazard(double x) const{
    return pow(x/scale_, shape_);
}

double hesim::Weibull::random() const{
    return R::rweibull(shape_, scale_);
}

/******************************
* Weibull distribution for NMA
******************************/
hesim::Weibull hesim::WeibullNma::create_from_Nma(double a0, double a1){
  double shape = a1 + 1;
  double scalePH = exp(a0)/shape;
  double scale = pow(scalePH, -1/shape);
  return Weibull(shape, scale);
}

hesim::WeibullNma::WeibullNma(double a0, double a1)
  : wei_(create_from_Nma(a0, a1)){
}

void hesim::WeibullNma::set_params(std::vector<double> params){
  wei_ = create_from_Nma(params[0], exp(params[1]));
}

double hesim::WeibullNma::pdf(double x) const{
    return wei_.pdf(x);
}

double hesim::WeibullNma::cdf(double x) const{
    return wei_.cdf(x);
}

double hesim::WeibullNma::quantile(double p) const{
    return wei_.quantile(p);
}

double hesim::WeibullNma::hazard(double x) const{
    return wei_.hazard(x);
}

double hesim::WeibullNma::cumhazard(double x) const{
    return wei_.cumhazard(x);
}

double hesim::WeibullNma::random() const{
    return wei_.random();
}

/*******************
* Gamma distribution
*******************/
hesim::Gamma::Gamma(double shape, double rate){
    shape_ = shape;
    rate_ = rate;
}

void hesim::Gamma::set_params(std::vector<double> params){
  shape_ = exp(params[0]);
  rate_ = exp(params[1]);
}

double hesim::Gamma::pdf(double x) const{
    return R::dgamma(x, shape_, 1/rate_, 0);
}

double hesim::Gamma::cdf(double x) const{
    return R::pgamma(x, shape_, 1/rate_, 1, 0);
}

double hesim::Gamma::quantile(double p) const{
    return R::qgamma(p, shape_, 1/rate_, 1, 0);
}

double hesim::Gamma::hazard(double x) const{
    return Gamma::pdf(x)/(1 - Gamma::cdf(x));
}

double hesim::Gamma::cumhazard(double x) const{
    return -R::pgamma(x, shape_, 1/rate_, 0, 1);
}

double hesim::Gamma::random() const{
    return R::rgamma(shape_, 1/rate_);
}

/************************
* Lognormal distribution
***********************/
hesim::Lognormal::Lognormal(double meanlog, double sdlog){
    meanlog_ = meanlog;
    sdlog_ = sdlog;
}

void hesim::Lognormal::set_params(std::vector<double> params){
  meanlog_ = params[0];
  sdlog_ = exp(params[1]);
}

double hesim::Lognormal::pdf(double x) const{
    return R::dlnorm(x, meanlog_, sdlog_, 0);
}

double hesim::Lognormal::cdf(double x) const{
    return R::plnorm(x, meanlog_, sdlog_, 1, 0);
}

double hesim::Lognormal::quantile(double p) const{
    return R::qlnorm(p, meanlog_, sdlog_, 1, 0);
}

double hesim::Lognormal::hazard(double x) const{
    return Lognormal::pdf(x)/(1 - Lognormal::cdf(x));
}

double hesim::Lognormal::cumhazard(double x) const{
    return -R::plnorm(x, meanlog_, sdlog_, 0, 1);
}

double hesim::Lognormal::random() const{
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

hesim::Gompertz::Gompertz(double shape, double rate){
    shape_ = shape;
    rate_ = rate;
}

void hesim::Gompertz::set_params(std::vector<double> params){
  shape_ = params[0];
  rate_ = exp(params[1]);
}

double hesim::Gompertz::pdf(double x) const{
    if (shape_ == 0){
        return R::dexp(x, 1/rate_, 0);
    }
    else{
        return rate_ * exp(shape_ * x) * exp(-rate_/shape_ * (exp(shape_ * x) -1));
    }
}

double hesim::Gompertz::cdf(double x) const{
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

double hesim::Gompertz::quantile(double p) const{
    return qgompertz(p, shape_, rate_);
}

double hesim::Gompertz::hazard(double x) const{
    return rate_ * exp(shape_ * x);
}

double hesim::Gompertz::cumhazard(double x) const{
    if (shape_ == 0){
        return rate_ * x;
    }
    else{
        return rate_/shape_ * expm1(shape_ * x);
    }
}

double hesim::Gompertz::random() const{
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

hesim::LogLogistic::LogLogistic(double shape, double scale){
    shape_ = shape;
    scale_ = scale;
}

void hesim::LogLogistic::set_params(std::vector<double> params){
  shape_ = exp(params[0]);
  scale_ = exp(params[1]);
}

double hesim::LogLogistic::pdf(double x) const{
    return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    pow((1 + pow((x/scale_), shape_)), 2);
}

double hesim::LogLogistic::cdf(double x) const{
    return 1 - 1/(1 + pow(x/scale_, shape_));
}

double hesim::LogLogistic::quantile(double p) const{
    return exp(R::qlogis(p, log(scale_), 1/shape_, 1, 0));
}

double hesim::LogLogistic::hazard(double x) const{
    return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    (1 + pow((x/scale_), shape_));
}

double hesim::LogLogistic::cumhazard(double x) const{
    return -log(1 - LogLogistic::cdf(x));
}

double hesim::LogLogistic::random() const{
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

hesim::GeneralizedGamma::GeneralizedGamma(double mu, double sigma, double Q){
    mu_ = mu;
    sigma_ = sigma;
    Q_ = Q;
}

void hesim::GeneralizedGamma::set_params(std::vector<double> params){
  mu_ = params[0];
  sigma_ = exp(params[1]);
  Q_ = exp(params[2]);
}

double hesim::GeneralizedGamma::pdf(double x) const{
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

double hesim::GeneralizedGamma::cdf(double x) const{
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

double hesim::GeneralizedGamma::quantile(double p) const{
    if (Q_ == 0){
        return R::qlnorm(p, mu_, 1/(sigma_ * sigma_), 1, 0);
    }
    else {
        double gamma_quantile = R::qgamma(p, 1/(Q_ * Q_), 1, 1, 0);
        return exp(mu_ + sigma_ * (log(Q_ * Q_ * gamma_quantile)/Q_));
    }
}

double hesim::GeneralizedGamma::hazard(double x) const{
    return GeneralizedGamma::pdf(x)/(1 - GeneralizedGamma::cdf(x));
}

double hesim::GeneralizedGamma::cumhazard(double x) const{
    return -log(1 - GeneralizedGamma::cdf(x));
}

double hesim::GeneralizedGamma::random() const{
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
        hesim::GeneralizedGamma gengamma(mu[index], sigma[index], Q[index]);
        sample[i] = gengamma.random();
    }
    return(sample);
}

/************************
* Royston/Parmar Splines
*************************/
hesim::SurvSplines::SurvSplines(std::vector<double> gamma,
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

void hesim::SurvSplines::set_params(std::vector<double> params){
  gamma_ = params;
}

double hesim::SurvSplines::timescale_fun(double x) const{
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

double hesim::SurvSplines::timescale_dx_fun(double x) const{
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

double hesim::SurvSplines::basis_cube(double x) const{
    if (x <= 0) {
        return 0;
    }
    else {
        return x * x * x;
    }
}

double hesim::SurvSplines::basis_cube_dx(double x) const{
    if (x <= 0) {
        return 0;
    }
    else {
        return 3 * x * x;
    }
}

double hesim::SurvSplines::linear_predict(double x) const{
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

double hesim::SurvSplines::linear_predict_dx(double x) const {
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

double hesim::SurvSplines::survival(double x) const{
    if (x <= 0){
        return 1; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard" | scale_ == "log_cumhazard"){
        return exp(-cumhazard(x));
    }
    else if (scale_ == "log_cumodds"){
        return 1/(1 + exp(linear_predict(x)));
    }
    else if (scale_ == "inv_normal"){
        return R::pnorm(-linear_predict(x), 0, 1, 1, 0);
    }
    else {
        Rcpp::stop("Selected scale is not available.");
    }
}

double hesim::SurvSplines::hazard(double x) const{
    if (x <= 0){
        return 0; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard"){
        return exp(linear_predict(x));
    }
    else if (scale_ == "log_cumhazard"){
        return timescale_dx_fun(x) * linear_predict_dx(x) * exp(linear_predict(x));
    }
    else if (scale_ == "log_cumodds"){
        return timescale_dx_fun(x) * linear_predict_dx(x) * R::plogis(linear_predict(x), 0, 1, 1, 0);
    }
    else if (scale_ == "inv_normal"){
        double lp = linear_predict(x);
        return timescale_dx_fun(x) * linear_predict_dx(x) * R::dnorm(-lp, 0, 1, 0)/R::pnorm(-lp, 0, 1, 1, 0);
    }
    else{
        Rcpp::stop("Selected scale is not available.");
    }
}

double hesim::SurvSplines::cumhazard(double x) const{
    if (x <= 0){
        return 0; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard"){
        return integrate_hazard(this, x);
    }
    else if (scale_ == "log_cumhazard"){
        return exp(linear_predict(x));
    }
    else if (scale_ == "log_cumodds"){
        return log1p(exp(linear_predict(x)));
    }
    else if (scale_ == "inv_normal"){
        return -R::pnorm(-linear_predict(x), 0, 1, 1, 1);
    }
    else{
        Rcpp::stop("Selected scale is not available.");
    }
}

double hesim::SurvSplines::pdf(double x) const {
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
    else if (scale_ == "log_cumodds"){
        prob = timescale_dx_fun(x) * linear_predict_dx(x) * exp(lp - 2 * log(1 + exp(lp)));
    }
    else if (scale_ == "inv_normal"){
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

double hesim::SurvSplines::cdf(double x) const{
    return 1 - survival(x);
}

double hesim::SurvSplines::quantile(double p) const{
  return quantile_numeric(this, p);
}

double hesim::SurvSplines::random() const{
    return quantile(R::runif(0, 1));
}

/************************
* Fractional polynomials
************************/
hesim::FracPoly::FracPoly(std::vector<double> gamma, std::vector<double> powers){
  gamma_ = gamma;
  powers_ = powers;
}

void hesim::FracPoly::set_params(std::vector<double> params){
  gamma_ = params;
}

double hesim::FracPoly::basis_power(double x, double power) const{
  if (power == 0){
    return log(x);
  }
  else{
    return pow(x, power);
  }
}

std::vector<double> hesim::FracPoly::basis(double x) const{
  int n_powers = powers_.size();
  std::vector<double> basis(n_powers + 1);
  basis[0] = 1;
  basis[1] = basis_power(x, powers_[0]);
  double xp_old = basis[1];
  double xp_new;
  if (n_powers > 1){
    for (int i = 1; i < n_powers; ++i){
      if (powers_[i] == powers_[i - 1]){
        xp_new = log(x) * xp_old;
      }
      else {
        xp_new = basis_power(x, powers_[i]);
      }
      basis[i + 1] = xp_new;
      xp_old = xp_new;
    }
  }
  return basis;
}

double hesim::FracPoly::linear_predict(double x) const{
  std::vector<double> b = basis(x);
  return std::inner_product(gamma_.begin(), gamma_.end(), b.begin(), 0.0);
}

double hesim::FracPoly::hazard(double x) const{
  if (x <= 0){
    return 0; //  model is for time >= 0
  }
  else{
    return exp(linear_predict(x));
  }
}

double hesim::FracPoly::cumhazard(double x) const{
  return integrate_hazard(this, x);
}

double hesim::FracPoly::cdf(double x) const{
  return 1 - exp(-cumhazard(x));
}

double hesim::FracPoly::pdf(double x) const{
  return hazard(x) * (1 - cdf(x));
}

double hesim::FracPoly::quantile(double p) const{
  return quantile_numeric(this, p);
}

double hesim::FracPoly::random() const{
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
    Rcpp::IntegerVector ans(k);
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


/**************
* RCPP Modules
**************/
RCPP_MODULE(Distributions){
  Rcpp::class_<hesim::Distribution>("Distribution")
  .method("pdf", &hesim::Distribution::pdf)
  .method("cdf", &hesim::Distribution::cdf)
  .method("quantile", &hesim::Distribution::quantile)
  .method("hazard", &hesim::Distribution::hazard)
  .method("cumhazard", &hesim::Distribution::cumhazard)
  .method("random", &hesim::Distribution::random)
  ;

  Rcpp::class_<hesim::Exponential>("Exponential")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double>()
    .method("pdf", &hesim::Exponential::pdf)
    .method("cdf", &hesim::Exponential::cdf)
    .method("quantile", &hesim::Exponential::quantile)
    .method("hazard", &hesim::Exponential::hazard)
    .method("cumhazard", &hesim::Exponential::cumhazard)
    .method("random", &hesim::Exponential::random)
  ;

  Rcpp::class_<hesim::Weibull>("Weibull")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::Weibull::pdf)
    .method("cdf", &hesim::Weibull::cdf)
    .method("quantile", &hesim::Weibull::quantile)
    .method("hazard", &hesim::Weibull::hazard)
    .method("cumhazard", &hesim::Weibull::cumhazard)
    .method("random", &hesim::Weibull::random)
  ;
  
  Rcpp::class_<hesim::WeibullNma>("WeibullNma")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::WeibullNma::pdf)
    .method("cdf", &hesim::WeibullNma::cdf)
    .method("quantile", &hesim::WeibullNma::quantile)
    .method("hazard", &hesim::WeibullNma::hazard)
    .method("cumhazard", &hesim::WeibullNma::cumhazard)
    .method("random", &hesim::WeibullNma::random)
  ;
  
  Rcpp::class_<hesim::Gamma>("Gamma")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::Gamma::pdf)
    .method("cdf", &hesim::Gamma::cdf)
    .method("quantile", &hesim::Gamma::quantile)
    .method("hazard", &hesim::Gamma::hazard)
    .method("cumhazard", &hesim::Gamma::cumhazard)
    .method("random", &hesim::Gamma::random)
  ;
  
  Rcpp::class_<hesim::Lognormal>("Lognormal")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::Lognormal::pdf)
    .method("cdf", &hesim::Lognormal::cdf)
    .method("quantile", &hesim::Lognormal::quantile)
    .method("hazard", &hesim::Lognormal::hazard)
    .method("cumhazard", &hesim::Lognormal::cumhazard)
    .method("random", &hesim::Lognormal::random)
  ;
  
  Rcpp::class_<hesim::Gompertz>("Gompertz")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::Gompertz::pdf)
    .method("cdf", &hesim::Gompertz::cdf)
    .method("quantile", &hesim::Gompertz::quantile)
    .method("hazard", &hesim::Gompertz::hazard)
    .method("cumhazard", &hesim::Gompertz::cumhazard)
    .method("random", &hesim::Gompertz::random)
  ;
  
  Rcpp::class_<hesim::LogLogistic>("LogLogistic")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::LogLogistic::pdf)
    .method("cdf", &hesim::LogLogistic::cdf)
    .method("quantile", &hesim::LogLogistic::quantile)
    .method("hazard", &hesim::LogLogistic::hazard)
    .method("cumhazard", &hesim::LogLogistic::cumhazard)
    .method("random", &hesim::LogLogistic::random)
  ;
  
  Rcpp::class_<hesim::GeneralizedGamma>("GeneralizedGamma")
    .derives<hesim::Distribution>("Distribution")
    .constructor<double, double, double>()
    .method("pdf", &hesim::GeneralizedGamma::pdf)
    .method("cdf", &hesim::GeneralizedGamma::cdf)
    .method("quantile", &hesim::GeneralizedGamma::quantile)
    .method("hazard", &hesim::GeneralizedGamma::hazard)
    .method("cumhazard", &hesim::GeneralizedGamma::cumhazard)
    .method("random", &hesim::GeneralizedGamma::random)
  ;
  
    Rcpp::class_<hesim::SurvSplines>("SurvSplines")
    .derives<hesim::Distribution>("Distribution")
    .constructor<std::vector<double>, std::vector<double>, std::string, std::string>()
    .method("linear_predict", &hesim::SurvSplines::linear_predict)
    .method("linear_predict_dx", &hesim::SurvSplines::linear_predict_dx)
    .method("pdf", &hesim::SurvSplines::pdf)
    .method("cdf", &hesim::SurvSplines::cdf)
    .method("quantile", &hesim::SurvSplines::quantile)
    .method("hazard", &hesim::SurvSplines::hazard)
    .method("cumhazard", &hesim::SurvSplines::cumhazard)
    .method("random", &hesim::SurvSplines::random)
  ;
    
     Rcpp::class_<hesim::FracPoly>("FracPoly")
    .derives<hesim::Distribution>("Distribution")
    .constructor<std::vector<double>, std::vector<double> >()
    .method("linear_predict", &hesim::FracPoly::linear_predict)
    .method("pdf", &hesim::FracPoly::pdf)
    .method("cdf", &hesim::FracPoly::cdf)
    .method("quantile", &hesim::FracPoly::quantile)
    .method("hazard", &hesim::FracPoly::hazard)
    .method("cumhazard", &hesim::FracPoly::cumhazard)
    .method("random", &hesim::FracPoly::random)
  ;
}
