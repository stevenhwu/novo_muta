/**
 * @file infinite_driver.cc
 * @author Melissa Ip
 *
 * This is a test driver file for using the expectation-maximization
 * implementation applied to the trio model.
 *
 * To compile on Herschel and include GSL:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o infinite_driver utility.cc read_dependent_data.cc trio_model.cc em_algorithm.cc sufficient_statistics.cc infinite_driver.cc
 */
#include "sufficient_statistics.h"


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
  SufficientStatistics stats(1.0);
  TrioVector sites = {data};
  stats.Update(params, sites);
  stats.Print();

  // M-Step.
  int count = 0;
  double maximized = stats.MaxSequencingErrorRate();
  while (!Equal(params.sequencing_error_rate(), maximized) && !stats.IsNan()) {  // Quits if converges.
    cout << "~E:\t" << params.sequencing_error_rate() << endl;
    cout << "~E':\t" << maximized << endl;
    params.set_sequencing_error_rate(maximized);  // Sets new estimate.
    stats.Clear(); // Sets to 0. 
    stats.Update(params, sites);  // Loops to E-Step.
    stats.Print();
    maximized = stats.MaxSequencingErrorRate(); // Loops to M-Step.
    count++;
  }
  
  cout << "^E:\t" << params.sequencing_error_rate() << endl;
  cout << count << " iterations." << endl;

  return 0;
}
