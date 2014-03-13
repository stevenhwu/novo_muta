/**
 * @file em_algorithm.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization
 * algorithm applied to a simplified and modified version of the trio model,
 * which will be updated as necessary to more complex and biologically realistic
 * models including the custom model using the Drichlet multinomial.
 */
#include "trio_model.cc"  // temp change to .h

// E-step methods.
double GetSomaticStatistic();
Matrix16_16d SomaticMutationCounts();


/**
 * Returns S_Som the number of nucleotide mismatches between all x and x′.
 * Calculated during the E-step of expectation-maximization algorithm.
 *
 * @return  number of expected somatic mutations
 */
double GetSomaticStatistic(TrioModel params) {
  // Matrix16_16d somatic_mutation_counts = SomaticMutationCounts();

  // RowVector16d s_som = RowVector16d::Zero();  // initially 0
  // RowVector16d s_som_mother = RowVector16d::Zero();
  // RowVector16d s_som_father = RowVector16d::Zero();
  // RowVector16d s_som_child = RowVector16d::Zero();

  // RowVector16d r_x = RowVector16d::Zero();
  // RowVector16d r_x_mother = RowVector16d::Zero();
  // RowVector16d r_x_father = RowVector16d::Zero();
  // RowVector16d r_x_child = RowVector16d::Zero();

  // double mother_term1 = 0.0;
  // double father_term1 = 0.0;
  // double child_term1 = 0.0;
  // double mother_term2 = 0.0;
  // double father_term2 = 0.0;
  // double child_term2 = 0.0;

  // // for each x genotype in population priors
  // for (int x = 0; x < kGenotypeCount * kGenotypeCount; ++x) {
  //   // P(R|somatic genotype)
  //   for (int i = 0; i < kGenotypeCount; ++i) {
  //     for (int j = 0; j < kGenotypeCount; ++j) {
  //       // s_som is 0, somatic_mutation_counts is term2
  //       mother_term1 = params.somatic_probability_mat()(i, j) * params.read_dependent_data()->mother_vec(j);
  //       father_term1 = params.somatic_probability_mat()(i, j) * params.read_dependent_data()->father_vec(j);
  //       child_term1 = params.somatic_probability_mat()(i, j) * params.read_dependent_data()->child_vec(j);

  //       r_x_mother(j) += mother_term1;
  //       r_x_father(j) += father_term1;
  //       r_x_child(j) += child_term1;

  //       // sum over y_j
  //       s_som_mother(j) += mother_term1 * somatic_mutation_counts(i, j) / mother_term1;
  //       s_som_father(j) += father_term1 * somatic_mutation_counts(i, j) / father_term1;
  //       s_som_child(j) += child_term1 * somatic_mutation_counts(i, j) / child_term1;
  //     }
  //   }

  //   // P(R|zygotic genotype)
  //   for (int i = 0; i < kGenotypeCount; ++i) {
  //     for (int j = 0; j < kGenotypeCount; ++j) {
  //       mother_term1 = (params.somatic_probability_mat()(i, j) *
  //         params.read_dependent_data()->denominator.mother_probability(j));
  //       father_term1 = (params.somatic_probability_mat()(i, j) *
  //         params.read_dependent_data()->denominator.father_probability(j));
  //       child_term1 = (params.somatic_probability_mat()(i, j) *
  //         params.read_dependent_data()->denominator.child_probability(j));

  //       r_x_mother(j) += mother_term1;
  //       r_x_father(j) += father_term1;
  //       r_x_child(j) += child_term1;

  //       mother_term2 = s_som_mother(j) + somatic_mutation_counts(i, j);
  //       father_term2 = s_som_father(j) + somatic_mutation_counts(i, j);
  //       child_term2 = s_som_child(j) + somatic_mutation_counts(i, j);

  //       // sum over y_j
  //       s_som_mother(j) += mother_term1 * mother_term2 / mother_term1; 
  //       s_som_father(j) += father_term1 * father_term2 / father_term1;
  //       s_som_child(j) += child_term1 * child_term2 / child_term1;
  //     }
  //   }

  //   // merge j branches
  //   s_som = s_som_mother + s_som_father;
  //   for (int i = 0; i < kGenotypeCount * kGenotypeCount; ++i) {
  //     for (int j = 0; j < kGenotypeCount; ++j) {
  //       child_term1 = (params.germline_probability_mat()(j, i) *
  //         params.read_dependent_data()->denominator.child_probability(j));
  //       s_som(j) += child_term1 * s_som_child(j) / child_term1;
  //     }
  //   }

  //   // calculate s_som at root of tree
  //   r_x = r_x_mother * r_x_father * r_x_child;
  // }

  return 0.0;  // temp
}

/**
 * Generates a 16 x 16 Eigen matrix holding the number of somatic mutations
 * where each allele site segregation is counted as a mutation in the infinite
 * sites model. For example:
 *
 *          Zygotic
 * Somatic  AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * AA       0  1  1  1  1  2  2  2  1  2  2  2  1  2  2  2
 * AC       1  0  1  1  2  1  2  2  2  1  2  2  2  1  2  2
 * ...            0
 *
 * @return  16 x 16 Eigen matrix holding the number of somatic mutations.
 */
Matrix16_16d SomaticMutationCounts() {
  Matrix16_16d mat = Matrix16_16d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      if (i != j) {  // leaves diagonal as 0's
        auto somatic_genotype = kGenotypeNumIndex.row(i);
        auto zygotic_genotype = kGenotypeNumIndex.row(j);
        // 1 mutation count per allele site segregation
        if (somatic_genotype(0) != zygotic_genotype(0)) {
          mat(i, j)++;
        }
        if (somatic_genotype(1) != zygotic_genotype(1)) {
          mat(i, j)++;
        }
      }
    }
  }
  return mat;
}