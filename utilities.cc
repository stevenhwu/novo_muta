/*
 * @file utilities.cc
 * @author Melissa Ip
 *
 * This contains useful constants and functions to support the TrioModel 
 * including the Dirichlet multinomial.
 */
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

// Global constants for specifying array size.
// Strings are lexicographically ordered.
const int kNucleotideCount = 4;
const char kNucleotides[kNucleotideCount] = {'A', 'C', 'G', 'T'};
const map<char, int> kNucleotideIndex = {
  {'A', 0},
  {'C', 1},
  {'G', 2},
  {'T', 3}
};
const int kGenotypeCount = 16;
const string kGenotypes[kGenotypeCount] = {
  "AA", "AC", "AG", "AT",
  "CA", "CC", "CG", "CT",
  "GA", "GC", "GG", "GT",
  "TA", "TC", "TG", "TT"
};
const map<string, int> kGenotypeIndex = {
  {"AA", 0},  {"AC", 1},  {"AG", 2},  {"AT", 3},
  {"CA", 4},  {"CC", 5},  {"CG", 6},  {"CT", 7},
  {"GA", 8},  {"GC", 9},  {"GG", 10}, {"GT", 11},
  {"TA", 12}, {"TC", 13}, {"TG", 14}, {"TT", 15}
};
const int kGenotypeNumIndex[kGenotypeCount][2] = {
  {0, 0}, {0, 1}, {0, 2}, {0, 3},
  {1, 0}, {1, 1}, {1, 2}, {1, 3},
  {2, 0}, {2, 1}, {2, 2}, {2, 3},
  {3, 0}, {3, 1}, {3, 2}, {3, 3}
};

union ReadData {
  uint16_t reads[kNucleotideCount];
  uint64_t key;
};

typedef Matrix<Array4d, kGenotypeCount, kGenotypeCount> Matrix4d1616;
typedef vector<ReadData> ReadDataVector;


/*
 * Initializes and returns a 16 x 16 x 4 Eigen matrix filled with zeros. The
 * third dimension is Array4d.
 *
 * @return 16 x 16 x 4 Eigen matrix filled with zeros.
 */
Matrix4d1616 ZeroMatrix4d1616() {
  Matrix4d1616 m;
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; i < kGenotypeCount; ++i) {
      m(i, j) = Array4d::Zero();
    }
  }
  return m;
}

/*
 * Prints a Matrix4d1616 object.
 *
 * @param m Matrix4d1616 object to be printed.
 */
void PrintMatrix4d1616(Matrix4d1616 m) {
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; i < kGenotypeCount; ++i) {
      cout << m(i, j) << endl;
    }
  }
}

/*
 * Generates a 16 x 4 alpha frequencies array given the sequencing error rate.
 * The order of the alpha frequencies is the same of that of kGenotypes.
 * 
 * Current values are placeholders until they are estimated in Spring 2014.
 *
 * @param rate Double representing sequencing error rate.
 * @return 16 x 4 Eigen array of Dirichlet multinomial alpha parameters
 *         alpha = (alpha_1, ..., alpha_K) for a K-category Dirichlet
 *         distribution (where K = 4 = kNucleotideCount) that vary with each
 *         combination of parental genotype and reference nucleotide.
 */
ArrayXXd GetAlphas(double rate) {
  ArrayXXd alphas(kGenotypeCount, kNucleotideCount);
  //        A             C             G             T
  alphas << 1 - rate,     rate/3,       rate/3,       rate/3,
            0.5 - rate/3, 0.5 - rate/3, rate/3,       rate/3,
            0.5 - rate/3, rate/3,       0.5 - rate/3, rate/3,
            0.5 - rate/3, rate/3,       rate/3,       0.5 - rate/3,

            0.5 - rate/3, 0.5 - rate/3, rate/3,       rate/3,
            rate/3,       1 - rate,     rate/3,       rate/3,
            rate/3,       0.5 - rate/3, 0.5 - rate/3, rate/3,
            rate/3,       0.5 - rate/3, rate/3,       0.5 - rate/3,

            0.5 - rate/3, rate/3,       0.5 - rate/3, rate/3,
            rate/3,       0.5 - rate/3, 0.5 - rate/3, rate/3,
            rate/3,       rate/3,       1 - rate,     rate/3,
            rate/3,       rate/3,       0.5 - rate/3, 0.5 - rate/3,

            0.5 - rate/3, rate/3,       rate/3,       0.5 - rate/3,
            rate/3,       0.5 - rate/3, rate/3,       0.5 - rate/3,
            rate/3,       rate/3,       0.5 - rate/3, 0.5 - rate/3,
            rate/3,       rate/3,       rate/3,       1 - rate;
  return alphas;
}

/*
 * Calculates probability from the probability density function (pdf):
 *
 * \frac{\gamma(\theta)}{\gamma(\theta + N)} *
 *     \Pi_{i = A, C, G, T} \frac{\gamma(\alpha_i * \theta + n_i)}
 *                               {\gamma(\alpha_i * \theta}
 *
 * @param alpha Eigen array of doubles representing frequencies for each
 *        category in the Dirichlet multinomial that sum to one.
 * @param n Array of integers representing samples for each category in the
 *          multinomial.
 * @return Double representing log_e(P) where P is the value calculated from the
 *         pdf.
 */   
double DirichletMultinomialLog(ArrayXd alpha, ReadData data) {
  double a = alpha.sum();
  int n = data.reads[0] + data.reads[1] + data.reads[2] + data.reads[3];
  double constant_term = lgamma(a) - lgamma(n + a);
  double product_term = 0.0;
  for (int i = 0; i < kNucleotideCount; ++i) {
    product_term += lgamma(alpha(i) + data.reads[i]) - lgamma(alpha(i));
  }
  return constant_term + product_term;
}

/*
 * Generates a 16 x 16 x 4 Eigen array. The first dimension represents the
 * mother genotype. The second dimension represents the father genotype. The
 * third dimension represents the nucleotide counts.
 *
 * Used in creating the population priors.
 *
 * @return 16 x 16 x 4 Eigen array of genotype counts where the (i, j) element
 *         is the count of nucleotide k of mother genotype i and father genotype
 *         j. 
 */
Matrix4d1616 TwoParentCounts() {
  Matrix4d1616 genotype_count = ZeroMatrix4d1616();
  for (const auto& nucleotide_pair : kNucleotideIndex) {
    char nucleotide = nucleotide_pair.first;
    int nucleotide_idx = nucleotide_pair.second;

    for (const auto& genotype_pair : kGenotypeIndex) {
      string genotype = genotype_pair.first;
      int genotype_idx = genotype_pair.second;

      for (auto base = genotype.begin(); base != genotype.end(); ++base) {
        if (*base == nucleotide) {
          auto row_slice = genotype_count.row(genotype_idx);
          for (int i = 0; i < row_slice.size(); ++i) {
            row_slice(i)(nucleotide_idx)++;
          }

          auto col_slice = genotype_count.col(genotype_idx);
          for (int i = 0; i < col_slice.size(); ++i) {
            col_slice(i)(nucleotide_idx)++;
          }
        }
      }
    }
  }
  return genotype_count;
}

int main(int argc, char *argv[]) {
  ArrayXd aa_alpha = GetAlphas(0.005).row(0);
  cout << aa_alpha << endl;
  ReadData data1 = {30, 0, 0, 0};
  cout << data1.reads << endl;
  double probability = DirichletMultinomialLog(aa_alpha, data1);
  cout << probability << endl;
  Matrix4d1616 m = TwoParentCounts();
  //PrintMatrix4d1616(m);
  return 0;
}

