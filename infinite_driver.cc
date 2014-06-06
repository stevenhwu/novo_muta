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
  	{40, 0, 0, 0},
  	{40, 0, 0, 0},
  	{40, 0, 0, 0}
  };

  double probability = params.MutationProbability(data);
  double somatic = GetSomaticStatistic(params);
  double germline = GetGermlineStatistic(params);
  SequencingErrorEstimates estimates;
  estimates.s_e = GetMismatchStatistic(params);
  estimates.s_hom = GetHomozygousStatistic(params);
  estimates.s_het = GetHeterozygousStatistic(params);
  double max = MaxSequencingErrorRate(estimates);

  cout << "P(Mut):\t" << probability << endl;
  cout << "S_Som:\t" << somatic << endl;
  cout << "S_Germ:\t" << germline << endl;
  cout << "S_E:\t" << estimates.s_e << endl;
  cout << "S_Hom:\t" << estimates.s_hom << endl;
  cout << "S_Het:\t" << estimates.s_het << endl;
  cout << "^E:\t" << max << endl;

  return 0;
}
