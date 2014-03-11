/**
 * @file infinite.cc
 * @author Melissa Ip
 *
 * This is a test driver file for using the expectation-maximization
 * implementation applied to the trio model. 
 */
#include "trio_model.h"


int main(int argc, const char *argv[]) {
  TrioModel params;
  ReadDataVector data = {
  	{0, 30, 0, 0},
  	{30, 0, 0, 0},
  	{30, 0, 0, 0}
  };

  cout << params.MutationProbability(data) << endl;
  cout << params.GetSomaticStatistic() << endl;

  return 0;
}
