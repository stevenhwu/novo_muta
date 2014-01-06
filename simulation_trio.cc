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

  ofstream fout;
  fout.open(argv[1]);
  TrioModel params;
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  for (auto trio : trio_vec) {
  	fout << params.MutationProbability(trio) << "\n";
  }
  fout.close();

  return 0;
}
