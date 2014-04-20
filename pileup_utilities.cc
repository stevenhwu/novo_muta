/**
 * @file pileup_utilities.cc
 * @author Melissa Ip
 *
 * This file contains methods to parse and align the the sequences in .pileup
 * files. Pileup files must be in the format described here:
 *
 * http://samtools.sourceforge.net/pileup.shtml
 *
 * This can create a TrioModel object using the parsed sequencing reads,
 * and write the probability of mutation to a text file.
 */
#include <fstream>
#include <iterator>
#include <sstream>

#include "trio_model.cc"

// Any greater probability than this number is printed.
const double kThreshold = 0.01;


/**
 * Removes initial sequences that are not necessary, which contain a N
 * reference (a site that is not matched to any base). Assume that there are
 * additional sequences where the reference is not N. The ifstream parameter is
 * changed. Returns the first line where N is not the reference.
 *
 * @param  f Input stream.
 * @return   First line where N is not the reference.
 */
string TrimHeader(ifstream &f) {
  string line;
  int sequence = 0;
  int position = 0;
  char ref_nucleotide;
  while(getline(f, line)) {
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    stringstream str(line);
    str >> sequence;
    str >> position;
    str >> ref_nucleotide;
    if (ref_nucleotide != 'N') {
      return line;
    }
  }
  return "";  // ERROR: This should not happen.
}

/**
 * Parses pileup data into ReadData. Periods and commas match the
 * reference nucleotide.
 *
 * @param  line Read from a pileup file representing a single site sequence.
 * @return      ReadData.
 */
ReadData GetReadData(const string &line) {
  int sequence = 0;
  int position = 0;
  char ref_nucleotide;
  int num_aligned_reads = 0;
  string bases;

  stringstream str(line);
  str >> sequence;
  str >> position;
  str >> ref_nucleotide;
  str >> num_aligned_reads;
  str >> bases;

  uint16_t matches = (count(bases.begin(), bases.end(), '.') +
                      count(bases.begin(), bases.end(), ','));
  uint16_t A = (count(bases.begin(), bases.end(), 'A') + 
                count(bases.begin(), bases.end(), 'a'));
  uint16_t C = (count(bases.begin(), bases.end(), 'C') + 
                count(bases.begin(), bases.end(), 'c'));
  uint16_t G = (count(bases.begin(), bases.end(), 'G') + 
                count(bases.begin(), bases.end(), 'g'));
  uint16_t T = (count(bases.begin(), bases.end(), 'T') + 
                count(bases.begin(), bases.end(), 't'));

  ReadData data;
  if (ref_nucleotide == 'A') {
    data = {matches, C, G, T};  
  } else if (ref_nucleotide == 'C') {
    data = {A, matches, G, T}; 
  } else if (ref_nucleotide == 'G') {
    data = {A, C, matches, T}; 
  } else if (ref_nucleotide == 'T') {
    data = {A, C, G, matches};
  }
  return data;
}

/**
 * Writes the probability of each site on a new line to a text file.
 *
 * @param  params      TrioModel object with default parameters.
 * @param  child_line  Line from the child pileup.
 * @param  mother_line Line from the mother pileup.
 * @param  father_line Line from the father pileup.
 */
double GetProbability(TrioModel &params, const string &child_line,
                      const string &mother_line, const string &father_line) {
  ReadDataVector data_vec = {GetReadData(child_line),
                             GetReadData(mother_line),
                             GetReadData(father_line)};
  return params.MutationProbability(data_vec);
}

/**
 * Opens and parses all pileup files. All valid sequences are converted to
 * ReadData and used to calculate the probability at their sequence position.
 * The output file is tab deliminated and each sequence is on a new line.
 * The first column represents the sequence position and the second column
 * represents the probability at that sequence.
 *
 * @param  file_name     Output file name.
 * @param  child_pileup  Chile pileup file name.
 * @param  mother_pileup Mother pileup file name.
 * @param  father_pileup Father pileup file name.
 */
void ProcessPileup(const string &file_name, const string &child_pileup,
                   const string &mother_pileup, const string &father_pileup) {
  ifstream child(child_pileup);
  ifstream mother(mother_pileup);
  ifstream father(father_pileup);
  if (!child.is_open() || 0 != child.fail() || !mother.is_open() ||
      0 != mother.fail() || !father.is_open() || 0 != father.fail()) {
    Die("Input file cannot be read.");
  }
  
  ofstream fout(file_name);
  TrioModel params;
  vector<double> probabilities;

  // Removes N sequences and writes probability of first valid line.
  string child_line = TrimHeader(child);
  string mother_line = TrimHeader(mother);
  string father_line = TrimHeader(father);
  if (child_line.empty() || mother_line.empty() || father_line.empty()) {
    Die("Pileup file does not contain valid sequences (no N reference).");
  }

  int sequence = 0;
  int position = 0;
  stringstream str(child_line);
  str >> sequence;
  str >> position;
  double probability = GetProbability(params, child_line, mother_line,
                                      father_line);
  if (probability >= kThreshold) {
    probabilities.push_back(probability);
  }

  // Writes probabilities of the rest of the sequences.
  while (getline(child, child_line)) {
    getline(mother, mother_line);
    getline(father, father_line);
    stringstream str(child_line);
    str >> sequence;
    str >> position;
    probability = GetProbability(params, child_line, mother_line, father_line);
    if (probability >= kThreshold) {
      probabilities.push_back(probability);
    }
  }

  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);

  child.close();
  mother.close();
  father.close();
  fout.close();
}
