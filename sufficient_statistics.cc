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
 *
 * @param  e              Expected statistic for sequencing error rate.
 * @param  hom            Expected statistic for homozygous matches.
 * @param  het            Expected statistic for heterozygous matches.
 * @param  som            Expected statistic for somatic mutation.
 * @param  germ           Expected statistic for germline mutation.
 * @param  log_likelihood Log likelihood of P(R,H)
 * @param  n_s            Number of sites.
 */
SufficientStatistics::SufficientStatistics(double sites_count)
    : e_{0.0}, hom_{0.0}, het_{0.0}, som_{0.0}, germ_{0.0}, log_likelihood_{0.0},
      n_s_{sites_count} {
}

/**
 * Maximizes germline mutation rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return  Maximized germline mutation rate.
 */
double SufficientStatistics::MaxGermlineMutationRate() {
  double bracket_term = 1.0 - (4.0/3.0) * ((2*germ_) / n_s_);
  return -0.75 * log(bracket_term);
}

/**
 * Maximizes somatic mutation rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return  Maximized somatic mutation rate.
 */
double SufficientStatistics::MaxSomaticMutationRate() {
  double bracket_term = 1.0 - (4.0/3.0) * ((6*som_) / n_s_);
  return -0.75 * log(bracket_term);
}

/**
 * Maximizes sequencing error rate. Calculated during the M-step of
 * expectation-maximization algorithm.
 *
 * @return  Maximized sequencing error rate.
 */
double SufficientStatistics::MaxSequencingErrorRate() {
  double sum = hom_ + het_ + e_;
  double double_het = 2.0 * het_;
  double sqrt_term_a = 9.0 * pow(hom_, 2.0);
  double sqrt_term_b = pow(double_het - e_, 2.0);
  double sqrt_term_c = 6.0 * hom_ * (double_het + e_);
  double sqrt_term = sqrt(sqrt_term_a + sqrt_term_b + sqrt_term_c);
  double inner_term = 3.0 * hom_ + double_het + 5.0 * e_;
  double subtract_term = (inner_term - sqrt_term) / sum / -3.0;
  double bracket_term = 1.0 + subtract_term;

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
    log_likelihood_ += log(params.read_dependent_data().denominator.sum);
  }
}

/**
 * Sets all statistics to 0 except n_s_.
 */
void SufficientStatistics::Clear() {
  e_ = 0.0;
  hom_ = 0.0;
  het_ = 0.0;
  som_ = 0.0;
  germ_ = 0.0;
  log_likelihood_ = 0.0;
}

/**
 * Returns true if any of the statistics is NaN.
 */
bool SufficientStatistics::IsNan() {
  return isnan(e_) || isnan(hom_) || isnan(het_) || isnan(som_) || isnan(germ_);
}

/**
 * Prints content.
 */
void SufficientStatistics::Print() {
  cout << "S_Som:\t"      << som_            << endl
       << "S_Germ:\t"     << germ_           << endl
       << "S_E:\t"        << e_              << endl
       << "S_Hom:\t"      << hom_            << endl
       << "S_Het:\t"      << het_            << endl
       << "Likelihood:\t" << log_likelihood_ << endl
       << "N_S:\t"        << n_s_            << endl;
}

double SufficientStatistics::e() const {
  return e_;
}

double SufficientStatistics::hom() const {
  return hom_;
}

double SufficientStatistics::het() const {
  return het_;
}

double SufficientStatistics::som() const {
  return som_;
}

double SufficientStatistics::germ() const {
  return germ_;
}

double SufficientStatistics::n_s() const {
  return n_s_;
}

double SufficientStatistics::log_likelihood() const {
  return log_likelihood_;
}

void SufficientStatistics::set_e(double max) {
  e_ = max;
}

void SufficientStatistics::set_hom(double max) {
  hom_ = max;
}

void SufficientStatistics::set_het(double max) {
  het_ = max;
}

void SufficientStatistics::set_som(double max) {
  som_ = max;
}

void SufficientStatistics::set_germ(double max) {
  germ_ = max;
}

void SufficientStatistics::set_n_s(double num) {
  n_s_ = num;
}

void SufficientStatistics::set_log_likelihood(double num) {
  log_likelihood_ = num;
}