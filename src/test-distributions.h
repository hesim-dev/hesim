# ifndef TEST_DISTRIBUTIONS_H
# define TEST_DISTRIBUTIONS_H

/**
 * \ingroup test
 * Templated function used to test hesim::stats::distribution.
 */
template <class Func>
double test_distribution(Func f, std::string fun, double x){
  if (fun == "pdf"){
     return f.pdf(x);
  }
  else if (fun == "cdf"){
    return f.cdf(x);
  }
  else if (fun == "quantile"){
    return f.quantile(x);
  }
  else if (fun == "hazard"){
    return f.hazard(x);
  }
  else if (fun == "cumhazard"){
    return f.cumhazard(x);
  }
  else if (fun == "random"){
    return f.random();
  }
  else {
    Rcpp::stop("Selected function is not available.");
  }  
}

# endif