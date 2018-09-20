// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <hesim/utils.h>
#include <hesim/stats/mean.h>

/**
 * Compute the incremental effect of treatment strategies (i.e., "interventions") 
 * relative to a comparator for an outcome of interest.
 * @param x Vector of outcomes for interventions.
 * @param y Vector of outcomes for the comparator. 
 * @param n_samples Number of random samples of the parameter values.
 * @param n_strategies Number of treatment strategies.
 * @param n_grps Number of subgroups.
 * @return A vector of incremental effects computed for each willingness to pay 
 * value in @p k, each of the @p n_samples randomly sampled parameter sets, 
 * and each of the @p n_grps subgroups. This value is exported to @c R and 
 * used in the internal function @c calc_incr_effect() (which is, in turn,
 * used in the functions @c incr_effect() and @c icea_pw()).
 */
// [[Rcpp::export]]
std::vector<double> C_incr_effect(std::vector<double> x, std::vector<double> y,
                                int n_samples, int n_strategies, int n_grps){
  int N = x.size();
  std::vector<double> incr_vec;
  incr_vec.reserve(N);

  // estimate incremental cost and qalys
  double incr = 0;
  int counter = 0;
  for (int g = 0; g < n_grps; ++g){
    for (int j = 0; j < n_strategies; ++j){
      for (int s = 0; s < n_samples; ++s){
        incr = x[counter] - y[g * n_samples + s];
        incr_vec.push_back(incr);
        ++counter;
      }
    }
  }
  return incr_vec;
}

// Probability CE in pairwise comparisons
/**
 * Compute values for cost-effectiveness acceptability curves (i.e., the 
 * probabilities that interventions are cost-effective relative to a comparator).
 * @param k Willingness to pay values.
 * @param ie Incremental clinical effectiveness.
 * @param ic Incremental costs.
 * @param n_samples Number of random samples of the parameter values.
 * @param n_strategies Number of treatment strategies.
 * @param n_grps Number of subgroups.
 * @return A vector where each element is the probability that the net 
 * monetary benefit of a given treatment strategy is greater than the net monetary
 * benefit of a comparator for a willingngess to pay value in @p k and each of 
 * the @p n_grps subgroups. This value is exported to @c R and used in the 
 * function @c ceac() (which is, in turn, used in @c icea()).
 */
// [[Rcpp::export]]
std::vector<double> C_ceac(std::vector<double> k, std::vector<double> ie,
                            std::vector<double> ic, int n_samples, int n_strategies, int n_grps) {

  // Initialize
  int n_k = k.size();
  int sumpos = 0;
  std::vector<double> prob;
  prob.reserve(n_k * n_strategies * n_grps);
  std::vector<double> k_vec;
  k_vec.reserve(n_k * n_strategies);
  double nb = 0;

  // loop over willingness to pay
  for (int j = 0; j < n_k; ++j){
    int counter = 0;

    // loop over groups
    for (int g = 0; g < n_grps; ++ g){

      // loop over treatment strategies
      for (int i = 0; i < n_strategies; ++i){
        sumpos = 0;

        // loop over simulations
        for (int s = 0; s < n_samples; ++s){
           nb = k[j] * ie[counter] - ic[counter];
          if (nb > 0.0){
            ++sumpos;
          }
          ++ counter;
        }
        prob.push_back((double)sumpos/n_samples);
        k_vec.push_back(k[j]);
      }
    }
  }
  return prob;
}

/**
 * Compute the probability that each treatment strategy is the most 
 * cost-effective (i.e., has the highest net monetary benefit).
 * @param k Willingness to pay values.
 * @param e Clinical effectiveness.
 * @param c Costs.
 * @param n_samples Number of random samples of the parameter values.
 * @param n_strategies Number of treatment strategies.
 * @param n_grps Number of subgroups.
 * @return A vector where each element is the probability that a given 
 * treatment strategy is the most cost-effective for a given willingness to pay 
 * value in @p k and one of the @p n_grps subgroups. This value is exported to 
 * @c R and used in the function @c mce() (which is, in turn, used in @c icea()).
 */
// [[Rcpp::export]]
std::vector<double> C_mce(std::vector<double> k, std::vector<double> e,
                 std::vector<double> c, int n_samples, int n_strategies, int n_grps) {

  int n_k = k.size();
  int N = n_k * n_strategies * n_grps;
  std::vector<double> prob(N, 0.0);
  int jg = 0;

  // loop over willingness to pay k
  for (int j = 0; j < n_k; ++j){
    int sg = 0;

    // loop over groups
    for (int g = 0; g < n_grps; ++g){

      // loop over simulations
      for (int s = 0; s < n_samples; ++s){
        std::vector<double> nb;
        nb.reserve(n_strategies);

        // loop over treatment strategies
        for (int i = 0; i < n_strategies; ++i){
          nb.push_back(k[j] * e[sg * n_strategies + i] - c[sg * n_strategies + i]);
        } // end strategies loop
        int nb_pos = hesim::max_element_pos(nb.begin(), nb.end());
        prob[jg * n_strategies + nb_pos] = prob[jg * n_strategies + nb_pos] + 1;
        ++sg;
      } // end samples loop
      ++jg;
    } // end group loop
  }

  // Convert sum to proportion
  for (int i = 0; i < N; ++i){
    prob[i] = prob[i]/n_samples;
  }

  return prob;
}

/**
 * Compute the expected net monetary benefit given perfect information (e.g.,
 * an average over the maximum net monetary benefit among available treatments
 *  at each simulation draw).
 * @param k Willingness to pay values.
 * @param e Clinical effectiveness.
 * @param c Costs.
 * @param n_samples Number of random samples of the parameter values.
 * @param n_strategies Number of treatment strategies.
 * @param n_grps Number of subgroups.
 * @return A vector of the expected net monetary benefit given perfect 
 * information computed for each willingness to pay value in @p k and each of the
 * @p n_grps subgroups. This value is exported to @c R and used in the 
 * function @c icea().
 */
// [[Rcpp::export]]
std::vector<double> C_enmbpi(std::vector<double> k,
                          std::vector<double> e, std::vector<double> c,
                          int n_samples, int n_strategies, int n_grps) {

  int n_k = k.size();
  std::vector<double> enmbpi;
  enmbpi.reserve(n_k * n_grps);

  // loop over willingness to pay k
  for (int j = 0; j < n_k; ++j){
    int sg = 0;

    // loop over groups
    for (int g = 0; g < n_grps; ++g){
      std::vector<double> nmbpi;
      nmbpi.reserve(n_samples);

      // loop over simulations
      for (int s = 0; s < n_samples; ++s){
        std::vector<double> inb;
        inb.reserve(n_strategies);

        // loop over treatment strategies
        for (int i = 0; i < n_strategies; ++i){
          inb.push_back(k[j] * e[sg * n_strategies + i] - c[sg * n_strategies + i]);
        }
        nmbpi.push_back(*std::max_element(inb.begin(), inb.end()));
        ++sg;
      }
      enmbpi.push_back(hesim::stats::mean(nmbpi));
    }
  }

  return enmbpi;
}

