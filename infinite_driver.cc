/**
 * @file infinite_driver.cc
 * @author Melissa Ip
 *
 * This is a test driver file for using the expectation-maximization
 * implementation applied to the trio model.
 *
 * To compile on Herschel and include GSL:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o infinite_driver utility.cc read_dependent_data.cc trio_model.cc em_algorithm.cc infinite_driver.cc
 */
#include "em_algorithm.h"


int main() {
  TrioModel params;
  ReadDataVector data = {
    {0, 0, 40, 0},
    {40, 0, 0, 0},
    {40, 0, 0, 0}
  };

  double probability = params.MutationProbability(data);
  cout << "P(Mut):\t" << probability << endl;

  // E-Step.
  ParamEstimates estimates;
  estimates.som = GetSomaticStatistic(params);
  estimates.germ = GetGermlineStatistic(params);
  estimates.e = GetMismatchStatistic(params);
  estimates.hom = GetHomozygousStatistic(params);
  estimates.het = GetHeterozygousStatistic(params);
  cout << "S_Som:\t"  << estimates.som  << endl
       << "S_Germ:\t" << estimates.germ << endl
       << "S_E:\t"    << estimates.e    << endl
       << "S_Hom:\t"  << estimates.hom  << endl
       << "S_Het:\t"  << estimates.het  << endl;
  
  // M-Step
  double maximized = MaxSequencingErrorRate(estimates);
  while (!Equal(params.sequencing_error_rate(), maximized)) {
    params.set_sequencing_error_rate(maximized);

    // Loops to E-Step.
    estimates.som = GetSomaticStatistic(params);
    estimates.germ = GetGermlineStatistic(params);
    estimates.e = GetMismatchStatistic(params);
    estimates.hom = GetHomozygousStatistic(params);
    estimates.het = GetHeterozygousStatistic(params);

    // Loops to M-Step. Quits if converges with sequencing error rate.
    maximized = MaxSequencingErrorRate(estimates);
  }

  cout << "^E:\t" << params.sequencing_error_rate() << endl;

  return 0;
}
