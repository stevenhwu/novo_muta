/**
 * @file em_algorithm.h
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization (EM)
 * algorithm.
 */
#ifndef EM_ALGORITHM_H
#define EM_ALGORITHM_H

#include "parameter_estimates.h"


/**
 * Returns ParameterEstimates object containing maximum likelihood estimates
 * and expected sufficient statistics for the sum of all sites by performing
 * the EM algorithm. Returns null if sites is empty.
 *
 * @param  params TrioModel object containing parameters.
 * @param  sites  List of trios.
 * @return        Parameter estimates.
 */
ParameterEstimates* EstimateParameters(TrioModel &params, const TrioVector &sites) {
  int sites_count = sites.size();
  if (sites_count > 0) {
    ParameterEstimates *stats = new ParameterEstimates(sites_count);
      // Exits if converges or takes longer than 50 iteratons.
      while (stats->Update(params, sites) &&
             !Equal(params.sequencing_error_rate(), stats->max_e()) &&
             stats->count() < 50) {
        params.set_sequencing_error_rate(stats->max_e());  // Sets new estimate.
      }
      return stats;
  }
  return NULL;
}

#endif
