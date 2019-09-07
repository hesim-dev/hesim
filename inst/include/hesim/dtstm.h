# ifndef HESIM_DTSTM_H
# define HESIM_DTSTM_H
#include <hesim/statmods/obs_index.h>
#include <hesim/statmods/statmods.h>

namespace hesim {

/** @ingroup dtstm 
 * Classes and functions for simulating discrete time state transiton models.
 */
namespace dtstm {

/***************************************************************************//** 
 * Simulate state probabilities
 * Simulate state probabilities using a cohort Markov model for a single patient.
 * @param x_first, x_last A container for storing the first and last values of
 * state probabilities for the patient. It is a stacked vector with each observation
 * denoting a health state and model cycle pair. 
 * @return 
 ******************************************************************************/ 
template <typename InputIt>
inline double sim_stateprobs(InputIt out_first, InputIt out_last){

}


} // end namespace ctstm

} // end namespace hesim

# endif
