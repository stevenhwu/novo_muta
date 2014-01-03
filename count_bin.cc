/**
 * @file count_bin.cc
 * @author Melissa Ip
 *
 * This script parses a simulation input file, where the first column represents
 * the probability of mutation (as a float [0, 1]) and the second column
 * represents whether the site contains a mutation (1 is true, 0 is false). Each
 * site is placed on a new line.
 *
 * This creates 10 bins numbered 0-9 with probability cateogories at 10%
 * intervals:
 * 
 * BIN   0        1         2        ...   9
 * %    [0, 10), [10, 20), [20, 30), ..., [90, 100]
 * 
 * This calculates the percentage of the sites in each bin that contain a real
 * mutation (value in the second column is 1). The digit in the tenths place of
 * the probability represents the number of the bin it belongs to. A probability
 * of 1.00 (100%) will go in the highest bin possible, bin 9.
 */
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

static const int kNumBins = 10;  // 10 bins cover 0-100% with 10% intervals.


/**
 * Prints user defined errors and exits the program.
 *
 * @param  msg Message to be printed.
 */
void Die(const char *msg) {
  if (errno == EDOM) {
    perror(msg);
  } else {
    printf("ERROR: %s\n", msg);
  }
  exit(EXIT_FAILURE);
}

int main(int argc, const char *argv[]) {
  if (argc < 2) {        
    Die("USAGE: count_bin <input.txt>");
  }

  const char *file_name = argv[1];
  string line;
  ifstream f(file_name);
  if (!f.is_open() || 0 != f.fail()) {
    Die("Input file cannot be read.");
  }

  double counts[kNumBins] = {0.0};
  double totals[kNumBins] = {0.0};
  double probability = 0.0;
  int has_mutation = 0;
  int bin = 0;
  int i = 0;

  while (getline(f, line)) {
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    stringstream str(line);
    str >> probability;
    str >> has_mutation;
    bin = (int) fmin(floor(probability * kNumBins), kNumBins - 1);
    totals[bin]++;
    if (has_mutation == 1) {
      counts[bin]++;
    }
  }
  f.close();

  for (i = 0; i < kNumBins; i++) {
    if (totals[i] > 0) {
      double has_mutation_percent = counts[i] / totals[i] * 100.0;
      printf("%.2f%% or %.f/%.f sites in bin %d contain a mutation.\n",
             has_mutation_percent, counts[i], totals[i], i);
    } else {
      printf("There are no sites in bin %d.\n", i);
    }
  }

  return 0;
}
