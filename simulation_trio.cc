/**
 * @file simulation_trio.cc
 * @author Melissa Ip
 *
 * This file outputs the probabilities of all possible unique trio sets at 4x
 * coverage.
 */
#include <fstream>

#include "trio_model.cc"


int main(int argc, const char *argv[]) {
  if (argc < 2) {
    Die("USAGE: simulation_trio <output.txt>");
  }

  const string file_name = argv[1];
  ofstream fout(file_name);
  TrioModel params;
  params.set_germline_mutation_rate(0.000001);
  params.set_somatic_mutation_rate(0.000001);
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  vector<double> probabilities;

  for (auto const &trio : trio_vec) {
    probabilities.push_back(params.MutationProbability(trio));
  }

  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);
  fout.close();

  return 0;
}
