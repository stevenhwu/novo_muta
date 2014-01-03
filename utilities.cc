/**
 * @file utilities.cc
 * @author Melissa Ip
 *
 * This file contains useful constants and functions to support the TrioModel 
 * including the Dirichlet multinomial and alphas.
 */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

union ReadData {
  uint16_t reads[4];
  uint64_t key;
};

typedef Array<int, 16, 2> Array16_2i;
typedef Array<double, 1, 16> Array16d;
typedef Array<double, 1, 256> Array256d;
typedef Array<double, 3, 16> Array3_16d;
typedef Array<double, 4, 16> Array4_16d;
typedef Array<double, 16, 4> Array16_4d;
typedef Array<double, 16, 16> Array16_16d;
typedef Array<double, 16, 256> Array16_256d;
typedef Matrix<Array4d, 16, 16> Matrix16_16_4d;
typedef vector<ReadData> ReadDataVector;

// Global constants for specifying array size and iterating through numeric
// representations of nucleotides and genotypes in lexicographical order.

// INDEX  GENOTYPE  NUCLEOTIDE 
// 0      AA        A
// 1      AC        C
// 2      AG        G
// 3      AT        T
// 4      CA
// 5      CC
// 6      CG
// 7      CT
// 8      GA
// 9      GC
// 10     GG
// 11     GT
// 12     TA
// 13     TC
// 14     TG
// 15     TT
static const int kGenotypeCount = 16;
static const int kNucleotideCount = 4;
static const double kEpsilon = numeric_limits<double>::epsilon();

/**
 * Returns 16 x 2 Eigen array where the first dimension represents genotypes
 * and the second dimension represents the numeric index of each nucleotide in
 * that genotype.
 *
 * Used for generating the constant kGenotypeNumIndex only. This constant is
 * used to determine what nucleotides a genotype is composed of or what
 * genotype a pair of nucleotides is.
 *
 * AA  {{0, 0},
 * AC   {0, 1},
 * AG   {0, 2},
 * AT   {0, 3},
 *       ...
 * TT   {3, 3}}
 *
 * @return  16 x 2 Eigen array.
 */
Array16_2i GenotypeNumIndex() {
  Array16_2i genotypeNumIndex;
  //                     A|    C|    G|    T|
  genotypeNumIndex << 0, 0, 0, 1, 0, 2, 0, 3,
                      1, 0, 1, 1, 1, 2, 1, 3,
                      2, 0, 2, 1, 2, 2, 2, 3,
                      3, 0, 3, 1, 3, 2, 3, 3;
  return genotypeNumIndex;
}
static const Array16_2i kGenotypeNumIndex = GenotypeNumIndex();


/**
 * Returns a 16 x 16 x 4 Eigen matrix filled with zeros. The third dimension is
 * Array4d.
 *
 * @return  16 x 16 x 4 Eigen matrix filled with zeros.
 */
Matrix16_16_4d ZeroMatrix16_16_4d() {
  Matrix16_16_4d m;
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; i < kGenotypeCount; ++i) {
      m(i, j) = Array4d::Zero();
    }
  }
  return m;
}

/**
 * Prints a Matrix16_16_4d object.
 *
 * @param  m Matrix16_16_4d object to be printed.
 */
void PrintMatrix16_16_4d(const Matrix16_16_4d& m) {
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; i < kGenotypeCount; ++i) {
      cout << m(i, j) << endl;
    }
  }
}

/**
 * Generates a 16 x 16 x 4 Eigen array. The first dimension represents the
 * mother genotype. The second dimension represents the father genotype. The
 * third dimension represents the nucleotide counts.
 *
 * Used for generating the constant kTwoParentCounts only, which is used in
 * creating the population priors.
 *
 * @return  16 x 16 x 4 Eigen array of genotype counts where the (i, j) element
 *          is the count of nucleotide k of mother genotype i and father
 *          genotype j.
 */
Matrix16_16_4d TwoParentCounts() {
  Matrix16_16_4d genotype_count = ZeroMatrix16_16_4d();
  // Loop through nucleotide and genotype indices.
  for (int nucleotide_idx = 0; nucleotide_idx < kNucleotideCount;
      ++nucleotide_idx) {
    for (int genotype_idx = 0; genotype_idx < kGenotypeCount; ++genotype_idx) {
      // Compare nucleotide with each nucleotide in genotype.
      for (int allele = 0; allele < 2; ++allele) {
        if (kGenotypeNumIndex(genotype_idx, allele) == nucleotide_idx) {
          // Increment nucleotide counts by 1 if match.
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
static const Matrix16_16_4d kTwoParentCounts = TwoParentCounts();

/**
 * Generates a 16 x 4 alpha frequencies array given the sequencing error rate.
 * The order of the alpha frequencies correspond to the genotypes. Each alpha
 * should sum to 1.
 * 
 * Current values are placeholders until they are estimated in Spring 2014.
 *
 * @param  rate Sequencing error rate.
 * @return      16 x 4 Eigen array of Dirichlet multinomial alpha parameters
 *              alpha = (alpha_1, ..., alpha_K) for a K-category Dirichlet
 *              distribution (where K = 4 = kNucleotideCount) that vary with
 *              each combination of parental genotype and reference nucleotide.
 */
Array16_4d GetAlphas(double rate) {
  Array16_4d alphas;
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

/**
 * Calculates probability from the probability density function (pdf):
 *
 * \frac{\gamma(\theta)}{\gamma(\theta + N)} *
 *     \Pi_{i = A, C, G, T} \frac{\gamma(\alpha_i * \theta + n_i)}
 *                               {\gamma(\alpha_i * \theta}
 *
 * @param  alpha Eigen array of doubles representing frequencies for each
 *               category in the Dirichlet multinomial.
 * @param  data  Read counts or samples for each category in the multinomial.
 * @return       log_e(P) where P is the value calculated from the pdf.
 */
double DirichletMultinomialLog(const Array4d& alpha, const ReadData& data) {
  double a = alpha.sum();
  int n = data.reads[0] + data.reads[1] + data.reads[2] + data.reads[3];
  double constant_term = lgamma(a) - lgamma(n + a);
  double product_term = 0.0;
  for (int i = 0; i < kNucleotideCount; ++i) {
    product_term += lgamma(alpha(i) + data.reads[i]) - lgamma(alpha(i));
  }
  return constant_term + product_term;
}

/**
 * Calculates the Kronecker product of an array with itself.
 *
 * Array sizes are specific to germline_probability_mat_.
 *
 * @param  arr 4 x 16 Eigen array.
 * @return     16 x 256 Eigen array.
 */
Array16_256d KroneckerProduct(const Array4_16d& arr) {
  Array16_256d kronecker_product;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 16; ++j) {
      for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 16; ++l) {
          kronecker_product(i*4 + k, j*16 + l) = arr(i, j) * arr(k, l);
        }
      }
    }
  }
  return kronecker_product;
}

/**
 * Calculates the Kronecker product of an array with itself.
 *
 * Array sizes are specific to somatic_probability_mat_.
 *
 * @param  arr 4 x 4 Eigen array.
 * @return     16 x 16 Eigen array.
 */
Array16_16d KroneckerProduct(const Array44d& arr) {
  Array16_16d kronecker_product;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 4; ++l) {
          kronecker_product(i*4 + k, j*4 + l) = arr(i, j) * arr(k, l);
        }
      }
    }
  }
  return kronecker_product;
}

/**
 * Calculates the Kronecker product of an array with itself.
 *
 * Array sizes are specific to the parent probability matrix.
 *
 * @param  arr1 1 x 16 Eigen array.
 * @param  arr2 1 x 16 Eigen array.
 * @return      16 x 256 Eigen array.
 */
Array256d KroneckerProduct(const Array16d& arr1, const Array16d& arr2) {
  Array256d kronecker_product;
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {
      kronecker_product(i*16 + j) = arr1(i) * arr2(j);
    }
  }
  return kronecker_product;
}

/**
 * Uses the dot product to multiply two arrays. Assume the arrays are the
 * appropriate dimensions.
 *
 * a x b arr1 *(dot) b x c arr2 = a x c product
 *
 * @param  arr1 Eigen array.
 * @param  arr2 Eigen array.
 * @return      Product of both arrays.
 */
ArrayXXd DotProduct(const ArrayXXd& arr1, const ArrayXXd& arr2) {
  int rows = arr1.rows();
  int cols = arr2.cols();
  ArrayXXd dot_product(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double sum = 0.0;
      for (int k = 0; k < arr1.cols(); ++k) {
        sum += arr1(i, k) * arr2(k, j);
      }
      dot_product(i, j) = sum;
    }
  }
  return dot_product;
}

/**
 * Returns 16 x 16 array with major diagonal unchanged and all other elements
 * replaced by 0. For example:
 *
 * {{1, 2, 3},      {{1, 0, 0},
 *  {4, 5, 6},  =>   {0, 5, 0},
 *  {7, 8, 9}}       {0, 0  9}}
 *
 * Array size is specific to somatic_probability_mat_.
 *
 * @param  arr 16 x 16 Eigen array.
 * @return     16 x 16 diagonal of array.
 */
Array16_16d GetDiagonal(const Array16_16d& arr) {
  Array16_16d diagonal = Array16_16d::Zero();
  auto vec = arr.matrix().diagonal();
  for (int i = 0; i < vec.size(); ++i) {
    diagonal(i, i) = vec(i);
  }
  return diagonal;
}

/**
 * Returns true if the two given doubles are equal to each other within epsilon
 * precision.
 *
 * @param  a First double.
 * @param  b Second double.
 * @return   True if the two doubles are equal to each other.
 */
bool Equal(double a, double b) {
  return fabs(a - b) < kEpsilon;
}

/**
 * Generates a random double in [min, max].
 *
 * @param  min Min range.
 * @param  max Max range.
 * @return     Random number in [min, max].
 */
double fRand(double min, double max) {
  double f = (double) rand() / RAND_MAX;
  return min + f * (max - min);
}

/**
 * Generates a random sample in the given the range.
 *
 * @param  range Samples generated from the set [0, range-1].
 * @param  p     Array of probabilities associated with each entry in the
 *               samples generated by range.
 * @return       Random element based on p probabilities from samples generated
 *               by range.
 */
int RandomChoice(int range, const ArrayXd &p) {
  // Picks random number using p probabilities.
  double weight_sum = p.sum();
  double random_idx = fRand(0.0, weight_sum);
  for (int i = 0; i < range; ++i) {
    if (random_idx < p(i)) {
      return i;
    }
    random_idx -= p(i);
  }
  return -1;  // ERROR: This should never happen.
}

/**
 * Generates a random sample in the given range.
 *
 * @param  range Samples generated from the set [0, range-1].
 * @param  p     Array of probabilities associated with each entry in the
 *               samples generated by range.
 * @param  size  Output shape, for generating multiple random samples.
 * @return       Random samples based on p probabilities from samples generated
 *               by range.
 */
ArrayXd RandomChoice(int range, const ArrayXd& p, int size) {
  // Creates size array to hold random samples.
  ArrayXd random_samples(size);
  while (size > 0) {
    size--;
    random_samples(size) = RandomChoice(range, p);
  }
  return random_samples;
}
