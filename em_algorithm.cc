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
  ReadDependentData *data = params.read_dependent_data();
  Matrix16_16d somatic_mutation_counts = SomaticMutationCounts();

  RowVector256d s_som = RowVector256d::Zero();  // initially 0
  RowVector16d s_som_mother = RowVector16d::Zero();
  RowVector16d s_som_father = RowVector16d::Zero();
  RowVector16d s_som_child = RowVector16d::Zero();
  RowVector256d s_som_mother_x = RowVector256d::Zero();
  RowVector256d s_som_father_x = RowVector256d::Zero();
  RowVector256d s_som_child_x = RowVector256d::Zero();

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;
  double mother_term2 = 0.0;
  double father_term2 = 0.0;
  double child_term2 = 0.0;

  // P(R|zygotic genotype)
  // S(R_mom, mom_zygotic=x), S(R_dad, dad_zygotic=x), S(R_child, child_zygotic=x)
  for (int x = 0; x < kGenotypeCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) {
      child_term1 = (params.somatic_probability_mat()(x, y) *
                     data->sequencing_probability_mat(0, y));
      mother_term1 = (params.somatic_probability_mat()(x, y) *
                      data->sequencing_probability_mat(1, y));
      father_term1 = (params.somatic_probability_mat()(x, y) *
                      data->sequencing_probability_mat(2, y));

      child_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      mother_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      father_term2 = /* 0 + */ somatic_mutation_counts(x, y);

      s_som_child(x) += child_term1 * child_term2;  // sum over y_j
      s_som_mother(x) += mother_term1 * mother_term2;
      s_som_father(x) += father_term1 * father_term2;
    }

    s_som_child(x) /= data->denominator.child_probability(x);
    s_som_mother(x) /= data->denominator.mother_probability(x);
    s_som_father(x) /= data->denominator.father_probability(x);
  }

  // for each genotype in population priors
  for (int x = 0; x < kGenotypeCount * kGenotypeCount; ++x) {
    // at top of branches
    // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
    // where s_som_mother and s_som_father do not change
    for (int y = 0; y < kGenotypeCount; ++y) {
      child_term1 = (params.germline_probability_mat()(y, x) *  // germline genotype given parent pair=x
                     data->denominator.child_probability(y));
      child_term2 = s_som_child(y) /* + 0 */;
      s_som_child_x(x) += child_term1 * child_term2;
    }

    s_som_child_x(x) /= data->denominator.child_germline_probability(x);
    
    // merge j branches
    // S(R_mom,R_dad,R_child, parent_pair=x)
    s_som(x) += (s_som_child_x(x) +
                 s_som_mother(x % kGenotypeCount) +
                 s_som_father(x / kGenotypeCount));

    // S(R_mom,R_dad,R_child)
    s_som(x) *= data->denominator.root_mat(x);
  }

  return s_som.sum() / data->denominator.sum;
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
