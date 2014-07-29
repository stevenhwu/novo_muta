/**
 * @file simulation_trio.cc
 * @author Melissa Ip
 *
 * This file outputs the probabilities of all possible trio sets at 4x
 * coverage. The population mutation rate is set to 0.001. The germline and
 * somatic mutation rates are both set to 1e-6.
 *
 * To compile on Herschel and include GSL:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o simulation_trio utility.cc read_dependent_data.cc trio_model.cc simulation_trio.cc
 *
 * To run this file, provide the following command line inputs:
 * ./simulation_trio <output.txt>
 */
#include <fstream>

#include "trio_model.h"


int main(int argc, const char *argv[]) {
  if (argc < 2) {
    Die("USAGE: simulation_trio <output.txt>");
  }

  const string file_name = argv[1];
  ofstream fout(file_name);
  TrioModel params;
  params.set_population_mutation_rate(0.001);
  params.set_germline_mutation_rate(1e-6);
  params.set_somatic_mutation_rate(1e-6);
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  vector<double> probabilities;

  for (const auto &trio : trio_vec) {
    probabilities.push_back(params.MutationProbability(trio));
  }

  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);
  fout.close();

  return 0;
}
