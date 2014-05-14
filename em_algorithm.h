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


// E-step methods.
RowVector16d GetMismatches(const ReadData &data);
double GetSequencingErrorStatistic(TrioModel params);
double GetGermlineStatistic(TrioModel params);
Matrix16_256d GermlineMutationCounts(TrioModel params);
Matrix4_16d GermlineMutationCountsSingle(TrioModel params);
double GetSomaticStatistic(TrioModel params);
Matrix16_16d SomaticMutationCounts();

#endif
