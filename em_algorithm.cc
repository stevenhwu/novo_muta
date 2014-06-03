/**
 * @file em_algorithm.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the expectation-maximization
 * algorithm.
 *
 * See top of em_algorithm.h for a complete description.
 */
#include "em_algorithm.h"


/**
 * Returns S_E the number of nucleotide mismatches between a somatic genotype
 * and its sequencing reads.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected sequencing errors.
 */
double GetSequencingErrorStatistic(const TrioModel& params) {
  ReadDependentData *data = params.read_dependent_data();
  const ReadDataVector &data_vec = data->read_data_vec;
  RowVector16d s_e_child = GetMismatches(data_vec[0]);
  RowVector16d s_e_mother = GetMismatches(data_vec[1]);
  RowVector16d s_e_father = GetMismatches(data_vec[2]);

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;

  // S(R_mom, mom_somatic=x), S(R_dad, dad_somatic=x), S(R_child, child_somatic=x)
  // for (int x = 0; x < kGenotypeCount; ++x) {
  //   for (int y = 0; y < kGenotypeCount; ++y) {
  //     child_term1 = data->child_somatic_probability(y);
  //     mother_term1 = data->mother_somatic_probability(y);
  //     father_term1 = data->father_somatic_probability(y);

  //     s_e_child(x) += child_term1;  // Sums over y_j.
  //     s_e_mother(x) += mother_term1;
  //     s_e_father(x) += father_term1;
  //   }
  // }

  return 0.0;
}

/**
 * Sums all nucleotide counts in ReadData and subtracts out the number of
 * nucleotides that match the genotype. It does not subtract twice for
 * homozygous genotypes. Returns 1 x 16 Eigen matrix holding number of
 * mismatches per genotype. For example:
 *
 * ReadData  A  C  G T
 *           20 10 0 1
 *
 * AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * 11 1  11 10 1  21 21 20 11 21 31 30 10 20 30 30
 *
 * @param  data ReadData.
 * @return      1 x 16 Eigen matrix containing number of mismatches per genotype.
 */
RowVector16d GetMismatches(const ReadData &data) {
  RowVector16d s_e = RowVector16d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    int allele1 = i / kNucleotideCount;
    int allele2 = i % kNucleotideCount;
    s_e(i) += data.reads[0] + data.reads[1] + data.reads[2] + data.reads[3];
    s_e(i) -= data.reads[allele1];  // Homozygous.
    if (allele1 != allele2) {
      s_e(i) -= data.reads[allele2];  // Hetereogyzous.
    }
  }

  return s_e;
}

/**
 * Returns 1 x 16 Eigen matrix holding number of heterozygous matches per
 * genotype. For example:
 *
 * ReadData  A  C  G T
 *           20 10 0 1
 *
 * AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * 0  30 20 21 30 0  10 11 20 10 0  1  21 11 1  0 
 *
 * @param  data ReadData.
 * @return      1 x 16 Eigen matrix containing number of heterozygous matches
 *              per genotype.
 */
RowVector16d GetHeterozygousMatches(const ReadData &data) {
  RowVector16d s_het = RowVector16d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    if (i % 5 != 0) {  // Homozygous genotypes are divisible by 5.
      s_e(i) += data.reads[i / kNucleotideCount] + data.reads[i % kNucleotideCount];
    }
  }

  return s_e;
}

/**
 * Returns 1 x 16 Eigen matrix holding number of homozygous matches per
 * genotype. For example:
 *
 * ReadData  A  C  G T
 *           20 10 0 1
 *
 * AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
 * 20 0  0  0  0  10 0  0  0  0  0  0  0  0  0  1 
 *
 * @param  data ReadData.
 * @return      1 x 16 Eigen matrix containing number of homozygous matches per
 *              genotype.
 */
RowVector16d GetHomozygousMatches(const ReadData &data) {
  RowVector16d s_hom = RowVector16d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    if (i % 5 == 0) {  // Homozygous genotypes are divisible by 5.
      s_e(i) += data.reads[i / kNucleotideCount];
    }
  }

  return s_e;
}

/**
 * Returns S_Germ the total number of nucleotide mismatches between parent and
 * child genotypes, or the expected number of germline mutations. This is the
 * sum of SM the number of nucleotide mismatches between m* and oa and SF the
 * number of nucleotide mismatches between f* and ob. Calculated during the
 * E-step of expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected germline mutations.
 */
double GetGermlineStatistic(const TrioModel& params) {
  ReadDependentData *data = params.read_dependent_data();
  Matrix16_256d germline_mutation_counts = GermlineMutationCounts(params);
  Matrix16_256d germline_probability_mat = params.germline_probability_mat();
  RowVector256d population_priors = params.population_priors();
  RowVector256d s_germ = RowVector256d::Zero();

  double child_term1 = 0.0;
  double child_term2 = 0.0;
  
  // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
  // At top of branches.
  // For each parent pair genotype x in population priors.
  for (int x = 0; x < kGenotypePairCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) {  // Child germline genotype.
      child_term1 = (germline_probability_mat(y, x) *
                     data->denominator.child_zygotic_probability(y));
      child_term2 = /* 0 + */ germline_mutation_counts(y, x);
      s_germ(x) += child_term1 * child_term2;
    }

    s_germ(x) *= data->denominator.mother_zygotic_probability(x / kGenotypeCount);
    s_germ(x) *= data->denominator.father_zygotic_probability(x % kGenotypeCount);
    s_germ(x) *= population_priors(x);
  }

  return s_germ.sum() / data->denominator.sum;
}

/**
 * Generates a 16 x 256 Eigen matrix holding the number of expected germline
 * mutations containing 0's, Het's, 1's, 1 + Het's, and 2's. The rows are child
 * germline genotypes and the columns are parent pair genotypes. Assumes that
 * the first genotype in the parent pair comes from the mother and the first
 * nucleotide of the child genotype comes from the mother, the second from the
 * father.
 *
 * @param  params TrioModel object containing parameters.
 * @return        16 x 256 Eigen matrix holding the number of germline mutations.
 */
Matrix16_256d GermlineMutationCounts(const TrioModel& params) {
  Matrix16_256d mat = Matrix16_256d::Zero();
  Matrix4_16d germline_mutation_counts = GermlineMutationCountsSingle(params);
  
  for (int y = 0; y < kGenotypeCount; ++y) {  // Child genotype.
    for (int x = 0; x < kGenotypePairCount; ++x) {  // Parent pair genotype.
      int child_allele1 = y / kNucleotideCount;  // Inherited from mother.
      int child_allele2 = y % kNucleotideCount;  // Inherited from father.
      int mother_genotype = x / kGenotypeCount;
      int father_genotype = x % kGenotypeCount;

      mat(y, x) = (germline_mutation_counts(child_allele1, mother_genotype) +
                   germline_mutation_counts(child_allele2, father_genotype));
    }
  }

  return mat;
}

/**
 * Generates a 4 x 16 Eigen matrix holding the number of expected germline
 * mutations or mismatches between the parent alleles and the zygote based on
 * a single parent. Assumes that the first genotype in the parent pair comes
 * from the mother and the first nucleotide of the child genotype comes from
 * the mother. For example:
 *
 *        Parent
 * Child  AA AC  AG  AT  CA  CC CG  CT  GA  GC  GG GT  TA  TC  TG  TT
 * A      0  Het Het Het Het 1  1   1   Het 1   1  1   Het 1   1   1
 * C      1  Het 1   1   Het 0  Het Het 1   Het 1  1   1   Het 1   1
 * G      1  1   Het 1   1   1  Het 1   Het Het 1  Het 1   1   Het 1
 * T      1  1   1   Het 1   1  1   Het 1   1   1  Het Het Het Het 1
 *
 * Het = 0.5 * no match / heterozygous match.
 *
 * @param  params TrioModel object containing parameters.
 * @return        4 x 16 Eigen matrix holding the number of germline mutations.
 */
Matrix4_16d GermlineMutationCountsSingle(const TrioModel& params) {
  Matrix4_16d mat = Matrix4_16d::Zero();

  for (int i = 0; i < kNucleotideCount; ++i) {  // Child allele.
    for (int j = 0; j < kGenotypeCount; ++j) {  // Single parent genotype.
      if (j % 5 == 0) {  // Homozygous parent genotypes are divisible by 5.
        if (IsAlleleInParentGenotype(i, j)) {
          mat(i, j) = 0.0;  // Allele is inherited, ex: A in AA.
        } else {
          mat(i, j) = 1.0;  // Allele is mutated, ex: C in AA.
        }
      } else {  // Heterozygous parent genotype.
        if (IsAlleleInParentGenotype(i, j)) {
          // Allele is potentially mutated (Het), ex: C in AC.
          mat(i, j) = 0.5 * params.no_match() / params.heterozygous_match();
        } else {
          mat(i, j) = 1.0;  // Allele is mutated, ex: A in CG.
        }
      }
    }
  }

  return mat;
}

/**
 * Returns S_Som the number of nucleotide mismatches between all x and x′,
 * or the expected number of somatic mutations. Calculated during the E-step of
 * expectation-maximization algorithm.
 *
 * @param  params TrioModel object containing parameters.
 * @return        Number of expected somatic mutations.
 */
double GetSomaticStatistic(const TrioModel& params) {
  ReadDependentData *data = params.read_dependent_data();
  Matrix16_16d somatic_probability_mat = params.somatic_probability_mat();
  Matrix16_16d somatic_mutation_counts = SomaticMutationCounts();
  Matrix16_256d germline_probability_mat = params.germline_probability_mat();

  RowVector256d s_som = RowVector256d::Zero();
  RowVector16d s_som_mother = RowVector16d::Zero();
  RowVector16d s_som_father = RowVector16d::Zero();
  RowVector16d s_som_child = RowVector16d::Zero();
  RowVector256d s_som_child_x = RowVector256d::Zero();  // Given x = parent pair genotype.

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;
  double mother_term2 = 0.0;
  double father_term2 = 0.0;
  double child_term2 = 0.0;

  // S(R_mom, mom_zygotic=x), S(R_dad, dad_zygotic=x), S(R_child, child_zygotic=x)
  // somatic_probability_mat and somatic_mutation_counts are both 16 x 16
  // matrices which are symmetrical across the major diagonal. Thus, this
  // assumes x represents zygotic while y represents somatic.
  for (int x = 0; x < kGenotypeCount; ++x) {  // Zygotic genotype.
    for (int y = 0; y < kGenotypeCount; ++y) {  // Somatic genotype.
      child_term1 = somatic_probability_mat(x, y) * data->child_somatic_probability(y);
      mother_term1 = somatic_probability_mat(x, y) * data->mother_somatic_probability(y);
      father_term1 = somatic_probability_mat(x, y) * data->father_somatic_probability(y);

      child_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      mother_term2 = /* 0 + */ somatic_mutation_counts(x, y);
      father_term2 = /* 0 + */ somatic_mutation_counts(x, y);

      s_som_child(x) += child_term1 * child_term2;  // Sums over y_j.
      s_som_mother(x) += mother_term1 * mother_term2;
      s_som_father(x) += father_term1 * father_term2;
    }

    // Dividing out child_zygotic_probability occurs at the next node.
    // s_som_child(x) /= data->denominator.child_zygotic_probability(x);
    s_som_mother(x) /= data->denominator.mother_zygotic_probability(x);
    s_som_father(x) /= data->denominator.father_zygotic_probability(x);
  }

  // S(R_mom, parent_pair=x), S(R_dad, parent_pair=x), S(R_child, parent_pair=x)
  // At top of branches where s_som_mother and s_som_father do not change,
  // because no somatic mutation can occur in germline.
  // For each parent pair genotype x in population priors.
  for (int x = 0; x < kGenotypePairCount; ++x) {
    for (int y = 0; y < kGenotypeCount; ++y) {  // Child germline genotype.
      child_term1 = germline_probability_mat(y, x);
      child_term2 = s_som_child(y) /* + 0 */;
      s_som_child_x(x) += child_term1 * child_term2;
    }

    s_som_child_x(x) /= data->denominator.child_germline_probability(x);  // Includes child_zygotic_probability.
    

    // S(R_mom,R_dad,R_child, parent_pair=x)
    // Merges j branches.
    s_som(x) += (s_som_child_x(x) +
                 s_som_mother(x / kGenotypeCount) +
                 s_som_father(x % kGenotypeCount));

    // S(R_mom,R_dad,R_child)
    // At root of tree.
    s_som(x) *= data->denominator.root_mat(x);  // Includes population_priors.
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
      if (i != j) {  // Diagonal is filled with 0's for no mutation.
        int somatic_allele1 = i / kNucleotideCount;
        int somatic_allele2 = i % kNucleotideCount;
        int zygotic_allele1 = j / kNucleotideCount;
        int zygotic_allele2 = j % kNucleotideCount;

        // 1 mutation count per allele site segregation.
        if (somatic_allele1 != zygotic_allele1) {
          mat(i, j)++;
        }
        if (somatic_allele2 != zygotic_allele2) {
          mat(i, j)++;
        }
      }
    }
  }

  return mat;
}
