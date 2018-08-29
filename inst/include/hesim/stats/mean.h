# ifndef HESIM_STATS_MEAN_H
# define HESIM_STATS_MEAN_H

# include <numeric>

namespace hesim {

namespace stats {

/**
 * @ingroup stats
 * Compute the arithmetic mean of a vector in
 * the Standard Library.
 * @param v A vector. Should be of type integer
 * or double.
 * @return None.
 */
inline double mean(std::vector<double> v) {
  return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

} // end namespace stats

} // end namespace hesim

# endif