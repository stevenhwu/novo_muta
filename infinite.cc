/**
 * @file infinite.cc
 * @author Melissa Ip
 *
 * This is a test driver file for using the expectation-maximization
 * implementation applied to the trio model.
 *
 * To compile on Herschel, use the following command to include the GSL library:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o infinite infinite.cc
 */
#include "em_algorithm.cc"


int main() {
  TrioModel params;
  // params.set_germline_mutation_rate(0.0000000002);  // smaller germline mutation rate increases likelihood of somatic mutation

  ReadDataVector data = {
  	{0, 0, 40, 0},
  	{40, 0, 0, 0},
  	{40, 0, 0, 0}
  };

  cout << params.MutationProbability(data) << endl;  // must call this or SetReadDependentData first to set ReadDependentData
  cout << GetSomaticStatistic(params) << endl;
  cout << GetGermlineStatistic(params) << endl;

  return 0;
}
