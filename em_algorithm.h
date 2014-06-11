/**
 * @file em_algorithm.h
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization
 * algorithm applied to a simplified and modified version of the trio model
 * (infinite sites model), which will be updated as necessary in Fall 2014 and
 * Spring 2015 to more complex and biologically realistic models including the
 * custom model using Dirichlet-multinomial approximations instead of
 * multinomial approximations.
 */
#ifndef EM_ALGORITHM_H
#define EM_ALGORITHM_H

#include "trio_model.h"

struct ParamEstimates {  // S_<T>
	double e;
	double hom;
	double het;
	double som;
	double germ;  // S_M + S_F.
	double n_s;  // Number of sites.
};

// M-step functions.
double MaxGermlineMutationRate(const ParamEstimates &estimates);
double MaxSomaticMutationRate(const ParamEstimates &estimates);
double MaxSequencingErrorRate(const ParamEstimates &estimates);

// E-step functions.
double GetPopulationMutationRateStatistic(const TrioModel &params);
double GetSequencingErrorStatistic(const TrioModel &params,
                                   const RowVector16d &child,
                                   const RowVector16d &mother,
                                   const RowVector16d &father);
double GetHomozygousStatistic(const TrioModel &params);
double GetHeterozygousStatistic(const TrioModel &params);
double GetMismatchStatistic(const TrioModel &params);
RowVector16d GetHeterozygousMatches(const ReadData &data);
RowVector16d GetHomozygousMatches(const ReadData &data);
RowVector16d GetMismatches(const ReadData &data);
double GetGermlineStatistic(const TrioModel &params);
Matrix16_256d GermlineMutationCounts(const TrioModel &params);
Matrix4_16d GermlineMutationCountsSingle(const TrioModel &params);
double GetSomaticStatistic(const TrioModel &params);
Matrix16_16d SomaticMutationCounts();

#endif
