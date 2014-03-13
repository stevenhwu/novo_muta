/**
 * @file infinite.cc
 * @author Melissa Ip
 *
 * This is a test driver file for using the expectation-maximization
 * implementation applied to the trio model. 
 */
#include "trio_model.cc"  // temp change to .h
//#include "em_algorithm.cc"

int main() {
  TrioModel params;
  ReadDataVector data = {
  	{0, 30, 0, 0},
  	{30, 0, 0, 0},
  	{30, 0, 0, 0}
  };

  cout << params.MutationProbability(data) << endl;
  // cout << GetSomaticStatistic(params) << endl;

  return 0;
}
