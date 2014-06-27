/**
 * @file parse_neg_trio.cc
 * @author Melissa Ip
 *
 * This outputs all of the trios that have negative probabilities from an input
 * file.
 */
#include <fstream>
#include <sstream>

#include "utility.h"


int main(int argc, const char *argv[]) {
  if (argc < 2) {        
    Die("USAGE: parse_neg_trio <input>.txt <output>.txt");
  }

  const string input = argv[1];
  const string output = argv[2];
  ifstream fin(input);
  if (!fin.is_open() || 0 != fin.fail()) {
    Die("Input file cannot be read.");
  }

  string line;
  int index = 0;
  double probability = 0.0;
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  
  ofstream fout(output);
  while (getline(fin, line)) {
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    stringstream str(line);
    str >> probability;

    if (probability < 0.0) {
      fout << trio_vec[index][0].reads[0] << " "
           << trio_vec[index][0].reads[1] << " "
           << trio_vec[index][0].reads[2] << " "
           << trio_vec[index][0].reads[3] << "\t"
           << trio_vec[index][1].reads[0] << " "
           << trio_vec[index][1].reads[1] << " "
           << trio_vec[index][1].reads[2] << " "
           << trio_vec[index][1].reads[3] << "\t"
           << trio_vec[index][2].reads[0] << " "
           << trio_vec[index][2].reads[1] << " "
           << trio_vec[index][2].reads[2] << " "
           << trio_vec[index][2].reads[3] << endl;
    }

    index++;
  }
  fin.close();
  fout.close();

  return 0;
}
