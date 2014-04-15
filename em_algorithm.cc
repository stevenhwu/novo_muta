/**
 * @file em_algorithm.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization
 * algorithm applied to a simplified and modified version of the trio model
 * (infinite sites model), which will be updated as necessary in Fall 2014 and
 * Spring 2015 to more complex and biologically realistic models including the
 * custom model using Dirichlet-multinomial approximations instead of
 * multinomial approximations.
 */
#include "trio_model.cc"  // FIXME: change to .h

// E-step methods.
double GetSomaticStatistic();
Matrix16_16d SomaticMutationCounts();


/**
 * Returns S_Som the number of nucleotide mismatches between all x and x′,
 * or the expected number of somatic mutations. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected somatic mutations.
 */
double GetSomaticStatistic(TrioModel params) {
  ReadDependentData *data = params.read_dependent_data();
  Matrix16_16d somatic_probability_mat = params.somatic_probability_mat();
  Matrix16_16d somatic_mutation_counts = SomaticMutationCounts();
  Matrix16_256d germline_probability_mat = params.germline_probability_mat();

  RowVector256d s_som = RowVector256d::Zero();
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
      child_term1 = somatic_probability_mat(x, y) * data->child_somatic_probability(y);
      mother_term1 = somatic_probability_mat(x, y) * data->mother_somatic_probability(y);
      father_term1 = somatic_probability_mat(x, y) * data->father_somatic_probability(y);

      child_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      mother_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      father_term2 = /* 0 + */ somatic_mutation_counts(x, y);

      s_som_child(x) += child_term1 * child_term2;  // sum over y_j
      s_som_mother(x) += mother_term1 * mother_term2;
      s_som_father(x) += father_term1 * father_term2;
    }

    s_som_child(x) /= data->denominator.child_zygotic_probability(x);
    s_som_mother(x) /= data->denominator.mother_zygotic_probability(x);
    s_som_father(x) /= data->denominator.father_zygotic_probability(x);
  }

  // for each genotype x in population priors
  for (int x = 0; x < kGenotypeCount * kGenotypeCount; ++x) {
    // at top of branches
    // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
    // where s_som_mother and s_som_father do not change
    for (int y = 0; y < kGenotypeCount; ++y) {
      child_term1 = (germline_probability_mat(y, x) *  // germline genotype given parent pair=x
                     data->denominator.child_zygotic_probability(y));  // uses zygotic
      child_term2 = s_som_child(y) /* + 0 */;
      s_som_child_x(x) += child_term1 * child_term2;
    }

    s_som_child_x(x) /= data->denominator.child_germline_probability(x);
    
    // merges j branches
    // S(R_mom,R_dad,R_child, parent_pair=x)
    s_som(x) += (s_som_child_x(x) + s_som_mother(x % kGenotypeCount) +
                 s_som_father(x / kGenotypeCount));

    // at root of tree
    // S(R_mom,R_dad,R_child)
    s_som(x) *= data->denominator.root_mat(x);  // root_mat includes population_priors
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
