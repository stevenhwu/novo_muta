/**
 * @file pileup_utility.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the functions needed to parse and
 * manipulate pileup data.
 *
 * See top of pileup_utility.h for a complete description.
 */
#include "pileup_utility.h"
 

/**
 * Removes initial sequences from the ifstream that are not necessary, which
 * contain a N reference (a site that is not matched to any base). Assumes that
 * all three pileup files contain the same number of invalid sequences. Returns
 * the first line where N is not the reference.
 *
 * @param  f Input stream.
 * @return   First line where N is not the reference.
 */
string PileupUtility::TrimHeader(ifstream &f) {
  string line;
  while(getline(f, line)) {
    line = GetSequence(line);  // Trims newline.
    if (!line.empty()) {
      return line;
    }
  }
  return "";  // ERROR: There were only invalid N sequences.
}

/**
 * Trims the newline fron the end of the line, and returns the line if it is a
 * valid sequence that does not contain a N reference.
 *
 * @param  line Line from pileup file.
 * @return      Line without newline and has valid nucleotide reference.
 */
string PileupUtility::GetSequence(string &line) {
  int sequence = 0;
  int position = 0;
  char ref_nucleotide;
  
  line.erase(remove(line.begin(), line.end(), '\n'), line.end());
  stringstream str(line);
  str >> sequence;
  str >> position;
  str >> ref_nucleotide;

  if (ref_nucleotide != 'N') {
    return line;
  } else {
    return "";
  }
}

/**
 * Parses pileup data into ReadData. Periods and commas match the
 * reference nucleotide.
 *
 * @param  line Read from a pileup file representing a single site sequence.
 * @return      ReadData.
 */
ReadData PileupUtility::GetReadData(const string &line) {
  int sequence = 0;
  int position = 0;
  char ref_nucleotide;
  int num_aligned_reads = 0;
  string bases;

  stringstream str(line);
  str >> sequence;
  str >> position;
  // str.seekg(10, ios::beg);
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
double PileupUtility::GetProbability(TrioModel &params, const string &child_line,
                      const string &mother_line, const string &father_line) {
  ReadDataVector data_vec = {GetReadData(child_line),
                             GetReadData(mother_line),
                             GetReadData(father_line)};
  return params.MutationProbability(data_vec);
}

/**
 * Opens and parses all pileup files. All valid sequences are converted to
 * ReadData and used to calculate the probability at their sequence position.
 * The output file is tab separated and each sequence is on a new line.
 * The first column represents the sequence position and the second column
 * represents the probability at that sequence.
 *
 * @param  file_name     Output file name.
 * @param  child_pileup  Chile pileup file name.
 * @param  mother_pileup Mother pileup file name.
 * @param  father_pileup Father pileup file name.
 */
void PileupUtility::WriteProbability(const string &file_name,
                                     const string &child_pileup,
                                     const string &mother_pileup,
                                     const string &father_pileup) {
  ifstream child(child_pileup);
  ifstream mother(mother_pileup);
  ifstream father(father_pileup);
  if (!child.is_open() || 0 != child.fail() || !mother.is_open() ||
      0 != mother.fail() || !father.is_open() || 0 != father.fail()) {
    Die("Input file cannot be read.");
  }
  
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

  child.close();
  mother.close();
  father.close();

  ofstream fout(file_name);
  ostream_iterator<double> output_iter(fout, "\n");
  copy(probabilities.begin(), probabilities.end(), output_iter);
  fout.close();
}
