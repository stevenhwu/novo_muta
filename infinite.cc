/**
 * @file infinite.cc
 * @author Melissa Ip
 *
 * This is a test driver file for using the expectation-maximization
 * implementation applied to the trio model.
 *
 * To compile on Herschel, use the following command to include the GSL library:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o infinite infinite.cc
 * 
 */
//#include "trio_model.h"  // temp change to .h
#include "em_algorithm.cc"

int main() {
  TrioModel params;
  ReadDataVector data = {
  	{20, 0, 20, 0},
  	{40, 0, 0, 0},
  	{40, 0, 0, 0}
  };

  cout << params.MutationProbability(data) << endl;
  cout << GetSomaticStatistic(params) << endl;

  return 0;
}
