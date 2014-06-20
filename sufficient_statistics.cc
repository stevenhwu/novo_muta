/**
 * @file sufficient_statistics.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the SufficientStatistics class.
 *
 * See top of sufficient_statistics.h for a complete description.
 */
#include "sufficient_statistics.h"


/**
 * Default constructor.
 */
SufficientStatistics::SufficientStatistics(double sites_count)
    : e_{0.0}, hom_{0.0}, het_{0.0}, som_{0.0}, germ_{0.0}, n_s_{sites_count} {
}

/**
 * Maximizes germline mutation rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return                  Maximized germline mutation rate.
 */
double SufficientStatistics::MaxGermlineMutationRate() {
  double bracket_term = 1.0 - 4.0/3.0 * germ_ / n_s_;
  return -0.75 * log(bracket_term);
}

/**
 * Maximizes somatic mutation rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return                  Maximized somatic mutation rate.
 */
double SufficientStatistics::MaxSomaticMutationRate() {
  double bracket_term = 1.0 - 4.0/3.0 * som_ / n_s_;
  return -0.75 * log(bracket_term);
}

/**
 * Maximizes sequencing error rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return                  Maximized sequencing error rate.
 */
double SufficientStatistics::MaxSequencingErrorRate() {
  double sum = hom_ + het_ + e_;
  double sqrt_term_a = 9.0 * pow(hom_, 2.0);
  double sqrt_term_b = pow(2.0 * het_ - e_, 2.0);
  double sqrt_term_c = 6.0 * hom_ * (2.0 * het_ + e_);
  double sqrt_term = sqrt(sqrt_term_a + sqrt_term_b + sqrt_term_c);
  double inner_term = 3.0 * hom_ + 2.0 * het_ + 5.0 * e_;
  double subtract_term = (inner_term - sqrt_term) / sum / 3.0;
  double bracket_term = 1.0 - subtract_term;
  return -0.75 * log(bracket_term);
}

/**
 * Calls the appropriate E-Step function for each summary statistic.
 *
 * @param params  TrioModel object containing parameters.
 * @param sites   List of parsed trios.
 */
void SufficientStatistics::Update(TrioModel &params, const TrioVector &sites) {
  for (const ReadDataVector data_vec : sites) {
    params.SetReadDependentData(data_vec);
    som_ += GetSomaticStatistic(params);
    germ_ += GetGermlineStatistic(params);
    e_ += GetMismatchStatistic(params);
    hom_ += GetHomozygousStatistic(params);
    het_ += GetHeterozygousStatistic(params);
  }
}

/**
 * Sets all statistics to 0 but keeps n_s_.
 */
void SufficientStatistics::Clear() {
  e_ = 0.0;
  hom_ = 0.0;
  het_ = 0.0;
  som_ = 0.0;
  germ_ = 0.0;
}

/**
 * Print content.
 */
void SufficientStatistics::Print() {
  cout << "S_Som:\t"  << som_  << endl
       << "S_Germ:\t" << germ_ << endl
       << "S_E:\t"    << e_    << endl
       << "S_Hom:\t"  << hom_  << endl
       << "S_Het:\t"  << het_  << endl
       << "N_S:\t"    << n_s_  << endl;
}
