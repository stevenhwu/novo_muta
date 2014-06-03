/*
 * @file counts_probability.cc
 * @author Melissa Ip
 *
 * This script parses a simulation input file, where the first column represents
 * the index of the trio in the reference TrioVector at 4x coverage, the
 * second column represents the number of random trios that have a mutation,
 * and the third column represents the number of random trios that do not have
 * a mutation. The rowwise sum of the second and third column is the total
 * number of random trios that match the key trio (given its index in the first
 * column).
 *
 * This calculates the empirical probability for each trio using the following
 * formula:
 *
 * P(mutation|trio) = #trios with mutation / #total trios
 *
 * The probabilities should match the probabilities from the result of
 * MutationProbability from simulation_trio.cc. This prints the probability
 * for each trio on a new line.
 */
#include <fstream>
#include <sstream>

#include "utility.h"


int main(int argc, const char *argv[]) {
  if (argc < 3) {
    Die("USAGE: counts_probability <input.txt> <output.txt>");
  }

  const char *file_name = argv[1];
  ifstream f(file_name);
  if (!f.is_open() || 0 != f.fail()) {
    Die("Input file cannot be read.");
  }

  const char *fout_name = argv[2];
  ofstream fout(fout_name);

  double probability = 0.0;
  vector<double> probabilities;
  int index = 0;  // Currently unused.
  int has_mutation_total = 0;
  int has_no_mutation_total = 0;
  int total_trios = 0;
  string line;

  while (getline(f, line)) {
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    stringstream str(line);
    str >> index;
    str >> has_mutation_total;
    str >> has_no_mutation_total;
    total_trios = has_mutation_total + has_no_mutation_total;
    if (total_trios == 0) {
      probability = 0.0;
    } else {
      probability = (double) has_mutation_total / total_trios;
    }
    probabilities.push_back(probability);
  }
  f.close();

  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);
  fout.close();

  return 0;
}
