/**
 * @file trio_model.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the TrioModel class.
 *
 * See top of trio_model.h for a complete description.
 */
#include "trio_model.h"


/**
 * Default constructor.
 * 
 * sequencing_probability_mat_ is created or updated if sequencing_error_rate_
 * or dirichlet_dispersion_ is changed when MutationProbability() is called.
 */
TrioModel::TrioModel()
    : population_mutation_rate_{0.001},
      germline_mutation_rate_{0.00000002},
      somatic_mutation_rate_{0.00000002},
      sequencing_error_rate_{0.005},
      dirichlet_dispersion_{1000.0},
      nucleotide_frequencies_{0.25, 0.25, 0.25, 0.25} {
  population_priors_ = TrioModel::PopulationPriors();
  germline_probability_mat_ = TrioModel::GermlineProbabilityMat();
  germline_probability_mat_num_ = TrioModel::GermlineProbabilityMat(true);
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMat();
  somatic_probability_mat_diag_ = TrioModel::SomaticProbabilityMatDiag();
  alphas_ = TrioModel::Alphas();
  read_dependent_data_.has_mutation = false;
}

/**
 * Constructor to customize parameters.
 *
 * @param  population_mutation_rate Population mutation rate.
 * @param  germline_mutation_rate   Germline mutation rate.
 * @param  somatic_mutation_rate    Somatic mutation rate.
 * @param  sequencing_error_rate    Sequencing error rate.
 * @param  dirichlet_dispersion     Dirichlet dispersion changes the magnitude
 *                                  of the probability of mutation.
 * @param  nucleotide_frequencies   Nucleotide frequencies changes the
 *                                  distribution of mutated nucleotides in the
 *                                  population priors.
 */
TrioModel::TrioModel(double population_mutation_rate,
                     double germline_mutation_rate,
                     double somatic_mutation_rate,
                     double sequencing_error_rate,
                     double dirichlet_dispersion,
                     const RowVector4d &nucleotide_frequencies)
    : population_mutation_rate_{population_mutation_rate},
      germline_mutation_rate_{germline_mutation_rate},
      somatic_mutation_rate_{somatic_mutation_rate},
      sequencing_error_rate_{sequencing_error_rate},
      dirichlet_dispersion_{dirichlet_dispersion},
      nucleotide_frequencies_{nucleotide_frequencies} {
  population_priors_ = TrioModel::PopulationPriors();
  germline_probability_mat_ = TrioModel::GermlineProbabilityMat();
  germline_probability_mat_num_ = TrioModel::GermlineProbabilityMat(true);
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMat();
  somatic_probability_mat_diag_ = TrioModel::SomaticProbabilityMatDiag();
  alphas_ = TrioModel::Alphas();
  read_dependent_data_.has_mutation = false;
}

/**
 * Returns S_Som the number of nucleotide mismatches between all x and x′.
 * Calculated during the E-step of expectation-maximization algorithm.
 *
 * @return  1 x 16 matrix containing number of expected somatic mutations
 */
 double TrioModel::GetSomaticStatistic() {
  Matrix16_16d somatic_mutation_counts = TrioModel::SomaticMutationCountsMatrix();

  RowVector16d s_som = RowVector16d::Zero();  // initially 0
  RowVector16d s_som_mother = RowVector16d::Zero();
  RowVector16d s_som_father = RowVector16d::Zero();
  RowVector16d s_som_child = RowVector16d::Zero();

  RowVector16d r_x = RowVector16d::Zero();
  RowVector16d r_x_mother = RowVector16d::Zero();
  RowVector16d r_x_father = RowVector16d::Zero();
  RowVector16d r_x_child = RowVector16d::Zero();

  double mother_term1 = 0.0;
  double father_term1 = 0.0;
  double child_term1 = 0.0;
  double mother_term2 = 0.0;
  double father_term2 = 0.0;
  double child_term2 = 0.0;

  // for each x genotype in population priors
  for (int x = 0; x < kGenotypeCount * kGenotypeCount; ++x) {
    // P(R|somatic genotype)
    for (int i = 0; i < kGenotypeCount; ++i) {
      for (int j = 0; j < kGenotypeCount; ++j) {
        // s_som is 0, somatic_mutation_counts is term2
        mother_term1 = somatic_probability_mat_(i, j) * read_dependent_data_.mother_vec(j);
        father_term1 = somatic_probability_mat_(i, j) * read_dependent_data_.father_vec(j);
        child_term1 = somatic_probability_mat_(i, j) * read_dependent_data_.child_vec(j);

        r_x_mother(j) += mother_term1;
        r_x_father(j) += father_term1;
        r_x_child(j) += child_term1;

        // sum over y_j
        s_som_mother(j) += mother_term1 * somatic_mutation_counts(i, j) / mother_term1;
        s_som_father(j) += father_term1 * somatic_mutation_counts(i, j) / father_term1;
        s_som_child(j) += child_term1 * somatic_mutation_counts(i, j) / child_term1;
      }
    }

    // P(R|zygotic genotype)
    for (int i = 0; i < kGenotypeCount; ++i) {
      for (int j = 0; j < kGenotypeCount; ++j) {
        mother_term1 = (somatic_probability_mat_(i, j) *
          read_dependent_data_.denominator.mother_probability(j));
        father_term1 = (somatic_probability_mat_(i, j) *
          read_dependent_data_.denominator.father_probability(j));
        child_term1 = (somatic_probability_mat_(i, j) *
          read_dependent_data_.denominator.child_probability(j));

        r_x_mother(j) += mother_term1;
        r_x_father(j) += father_term1;
        r_x_child(j) += child_term1;

        mother_term2 = s_som_mother(j) + somatic_mutation_counts(i, j);
        father_term2 = s_som_father(j) + somatic_mutation_counts(i, j);
        child_term2 = s_som_child(j) + somatic_mutation_counts(i, j);

        // sum over y_j
        s_som_mother(j) += mother_term1 * mother_term2 / mother_term1; 
        s_som_father(j) += father_term1 * father_term2 / father_term1;
        s_som_child(j) += child_term1 * child_term2 / child_term1;
      }
    }

    // merge j branches
    s_som = s_som_mother + s_som_father;
    for (int i = 0; i < kGenotypeCount * kGenotypeCount; ++i) {
      for (int j = 0; j < kGenotypeCount; ++j) {
        child_term1 = (germline_probability_mat_(j, i) *
          read_dependent_data_.denominator.child_probability(j));
        s_som(j) += child_term1 * s_som_child(j) / child_term1;
      }
    }

    // calculate s_som at root of tree
    r_x = r_x_mother * r_x_father * r_x_child;
  }

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
Matrix16_16d TrioModel::SomaticMutationCountsMatrix() {
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

/**
 * Implements the trio model for a single site by calling the functions on the
 * left in the following diagram. The function names label the arrow-denoted
 * processes in the population model.
 *
 *                             Population          Population
 * PopulationPriors()              |                   |
 *                                 v                   v
 *                               Mother              Father
 *                               Zygotic             Zygotic
 *                               Diploid             Diploid
 *                               Genotype            Genotype
 * GermlineMutation()              |   \             / |  
 * GermlineProbabilityMat()        |    v           v  |
 *                                 |      Daughter     |
 *                                 |      Diploid      |
 * SomaticMutation()               |         |         |
 * SomaticProbabilityMat()         v         v         v
 *                               Mother   Daughter  Father
 *                               Somatic  Somatic   Somatic
 *                               Diploid  Diploid   Diploid
 *                               Genotype Genotype  Genotype
 * SequencingProbabilityMat()      |         |         |
 *                                 v         v         v
 *                               Mother   Daughter  Father
 *                               Genotype Genotype  Genotype
 *                               Reads    Reads     Reads
 *
 * Initializes and updates read_dependent_data.sequencing_probability_mat,
 * read_dependent_data.child_vec, read_dependent_data.mother_vec, and
 * read_dependent_data.father_vec using sequencing_error_rate_.
 * dirichlet_dispersion_ is not used.
 *
 * @param   data_vec Read counts in order of child, mother and father.
 * @return           Probability of mutation given read data and parameters.
 */
double TrioModel::MutationProbability(const ReadDataVector &data_vec) {
  TrioModel::SequencingProbabilityMat(data_vec);
  TrioModel::SomaticTransition();
  TrioModel::GermlineTransition();
  TrioModel::SomaticTransition(true);
  TrioModel::GermlineTransition(true);

  return 1 - (read_dependent_data_.numerator.sum /
    read_dependent_data_.denominator.sum);
}

/**
 * Returns 1 x 256 Eigen probability RowVector. This is an order-relevant
 * representation of the possible events in the sample space that covers all
 * possible parent genotype combinations. For example:
 *
 * [P(AAAA), P(AAAC), P(AAAG), P(AAAT), P(AACA), P(AACC), P(AACG)...]
 *
 * Sets genotype_mat_ as the 1 x 16 Eigen probability RowVector for a single
 * parent.
 *
 * Calls SpectrumProbability assuming infinite sites model using all enumerated
 * nucleotide counts at coverage 4x.
 *
 * @return  1 x 256 Eigen probability RowVector in log e space where the (i, j)
 *          element is the probability that the mother has genotype i and the
 *          father has genotype j.
 */
RowVector256d TrioModel::PopulationPriors() {
  Matrix16_16d population_priors = Matrix16_16d::Zero();
  // Resizes population_priors to 1 x 256 RowVector.
  RowVector256d population_priors_flattened;
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      RowVector4d nucleotide_counts = kTwoParentCounts(i, j);
      double probability = TrioModel::SpectrumProbability(nucleotide_counts);
      population_priors(i, j) = probability;
      int idx = i * kGenotypeCount + j;
      population_priors_flattened(idx) = population_priors(i, j);
    }
  }

  // Sets genotype_mat_ to collapsed single parent probability matrix.
  genotype_mat_ = population_priors.rowwise().sum();

  return population_priors_flattened;
}

/**
 * Returns the probability of the drawn alleles having a 4-0, 3-1, or 2-2
 * spectrum.
 *
 * @param  nucleotide_counts RowVector containing allele counts.
 * @return                   Probability of allele spectrum.
 */
double TrioModel::SpectrumProbability(const RowVector4d &nucleotide_counts) {
  double p2_2 = population_mutation_rate_ / 2.0;
  double p3_1 = population_mutation_rate_ + population_mutation_rate_ / 3.0;
  double p4_0 = 1.0 - p3_1 - p2_2;

  if (IsInVector(nucleotide_counts, 4.0)) {
    return p4_0 * 0.25;  // p(4 allele) * p(position)
  } else if (IsInVector(nucleotide_counts, 3.0)) {
    return p3_1 * 0.25 / 3.0 * 0.25;  // p(3 allele) * p(1 allele) * p(position)
  } else if (IsInVector(nucleotide_counts, 2.0) &&
      !IsInVector(nucleotide_counts, 1.0)) {
    return p2_2 * 0.25 / 3.0 * 2.0 / 6.0;  // p(2 allele) * p(2 allele) * pair qualifier * p(position)
  } else {
    return 0.0;  // counts do not match 4-0, 3-0, or 2-2 allele spectrum
  }
}

/**
 * Calculates the probability of germline mutation and parent chromosome
 * donation by default. Assume the first chromosome is associated with
 * the mother and the second chromosome is associated with the father.
 *
 * @param  child_nucleotide_idx Index of child allele.
 * @param  parent_genotype_idx  Index of parent genotype.
 * @param  no_mutation_flag     False by default. Set to true to calculate
 *                              probability of no mutation.
 * @return                      Probability of germline mutation.
 */
double TrioModel::GermlineMutation(int child_nucleotide_idx,
                                   int parent_genotype_idx,
                                   bool no_mutation_flag) {
  // Calculates set of possible probabilities using mutation rate.
  double exp_term = exp(-4.0/3.0 * germline_mutation_rate_);
  double homozygous_match = 0.25 + 0.75 * exp_term;
  double heterozygous_match = 0.25 + 0.25 * exp_term;
  double no_match = 0.25 - 0.25 * exp_term;
  
  // Checks if the child nucleotide is in the parent genotype.
  auto parent_genotype = kGenotypeNumIndex.row(parent_genotype_idx);
  bool is_in_vec = false;
  if (child_nucleotide_idx == parent_genotype(0) ||
      child_nucleotide_idx == parent_genotype(1)) {
    is_in_vec = true;
  }

  // Determines if the comparison is homozygous, heterozygous or no match.
  if (is_in_vec) {
    if (parent_genotype(0) == parent_genotype(1)) {
      return homozygous_match;
    } else {
      if (no_mutation_flag) {
        return homozygous_match / 2;
      } else {
        return heterozygous_match;
      }
    }
  } else {
    if (no_mutation_flag) {
      return 0.0;
    } else {
      return no_match;
    }
  }
}

/**
 * Calculates the probability matrix for the offspring using the given mutation
 * rate, derived from the Kronecker product of a 4 x 16 Eigen parent 
 * probability matrix with itself.
 *
 * This is a transition matrix used to mutate the child germline genotype.
 *
 * @param  no_mutation_flag False by default. Set to true to calculate
 *                          probability of no mutation.
 * @return                  16 x 256 Eigen probability matrix.
 */
Matrix16_256d TrioModel::GermlineProbabilityMat(bool no_mutation_flag) {
  Matrix4_16d germline_probability_mat = Matrix4_16d::Zero();
  for (int i = 0; i < kNucleotideCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      double probability = TrioModel::GermlineMutation(i, j, no_mutation_flag);
      germline_probability_mat(i, j) = probability;
    }
  }
  return KroneckerProduct(germline_probability_mat);
}

/**
 * Calculates the probability of somatic mutation.
 *
 * @param  nucleotide_idx       Index of nucleotide.
 * @param  other_nucleotide_idx Index of another nucleotide to be compared.
 * @return                      Probability of somatic mutation.
 */
double TrioModel::SomaticMutation(int nucleotide_idx, int other_nucleotide_idx) {
  double exp_term = exp(-4.0/3.0 * somatic_mutation_rate_);
  double term = 0.25 - 0.25 * exp_term;
  double indicator_term = 0.0;

  // Checks if indicator function is true for each chromosome.
  if (nucleotide_idx == other_nucleotide_idx) {
    indicator_term = exp_term;
  }
  return term + indicator_term;
}

/**
 * Computes event space for somatic nucleotide given a genotype nucleotide for
 * a single chromosome. Combines event spaces for two chromosomes independent
 * of each other and calculates somatic mutation probability matrix for a 
 * single parent.
 *
 * This is a transition matrix used to mutate the somatic genotypes.
 *
 * @return  16 x 16 Eigen probability matrix where the first dimension is the
 *          original somatic genotypes and the second dimension is the mutated 
 *          genotype.
 */
Matrix16_16d TrioModel::SomaticProbabilityMat() {
  Matrix4d somatic_probability_mat = Matrix4d::Zero();
  for (int i = 0; i < kNucleotideCount; ++i) {
    for (int j = 0; j < kNucleotideCount; ++j) {
      double probability = TrioModel::SomaticMutation(i, j);
      somatic_probability_mat(i, j) = probability;
    }
  }
  return KroneckerProduct(somatic_probability_mat);
}

/**
 * Returns diagonal of somatic_probability_mat_. Assumes SomaticProbabilityMat
 * has been called.
 *
 * @return  16 x 16 Eigen matrix diagonal of somatic_probability_mat_.
 */
Matrix16_16d TrioModel::SomaticProbabilityMatDiag() {
  return somatic_probability_mat_.diagonal().asDiagonal();
}

/**
 * Calculates the probability of sequencing error for all read data. Assume
 * data contains 3 reads (child, mother, father). Assume each chromosome is
 * equally likely to be sequenced.
 *
 * Adds the max element of all reads in ReadDataVector to
 * read_dependent_data_.max_elements before rescaling to normal space.
 *
 * @param  data_vec ReadDataVector containing nucleotide counts for trio family.
 */
void TrioModel::SequencingProbabilityMat(const ReadDataVector &data_vec) {
  read_dependent_data_.sequencing_probability_mat = Matrix3_16d::Zero();
  for (int read = 0; read < 3; ++read) {
    for (int genotype_idx = 0; genotype_idx < kGenotypeCount; ++genotype_idx) {
      auto alpha = alphas_.row(genotype_idx);
      // converts alpha to double array
      double p[kNucleotideCount] = {alpha(0), alpha(1), alpha(2), alpha(3)};
      // converts read to unsigned int array
      ReadData data = data_vec[read];
      unsigned int n[kNucleotideCount] = {data.reads[0], data.reads[1],
                                          data.reads[2], data.reads[3]};
      double log_probability = gsl_ran_multinomial_lnpdf(kNucleotideCount, p, n);
      read_dependent_data_.sequencing_probability_mat(read, genotype_idx) = log_probability;
    }
  }
  
  // Rescales to normal space and records max element of all 3 reads together.
  double max_element = read_dependent_data_.sequencing_probability_mat.maxCoeff();
  read_dependent_data_.max_elements.push_back(max_element);

  // Calculates sequencing_probability_mat and splits into individual child
  // mother, and father vectors.
  read_dependent_data_.sequencing_probability_mat = exp(
    read_dependent_data_.sequencing_probability_mat.array() - max_element
  );
  read_dependent_data_.child_vec = read_dependent_data_.sequencing_probability_mat.row(0);
  read_dependent_data_.mother_vec = read_dependent_data_.sequencing_probability_mat.row(1);
  read_dependent_data_.father_vec = read_dependent_data_.sequencing_probability_mat.row(2);
}

/**
 * Multiplies sequencing probability vectors by somatic transition matrix.
 *
 * @param  is_numerator True if calculating probability of numerator.
 */
void TrioModel::SomaticTransition(bool is_numerator) {
  if (!is_numerator) {
    read_dependent_data_.denominator.child_probability = (
      read_dependent_data_.child_vec * somatic_probability_mat_
    );
    read_dependent_data_.denominator.mother_probability = (
      read_dependent_data_.mother_vec * somatic_probability_mat_
    );
    read_dependent_data_.denominator.father_probability = (
      read_dependent_data_.father_vec * somatic_probability_mat_
    );
  } else {
    read_dependent_data_.numerator.child_probability = (
      read_dependent_data_.child_vec * somatic_probability_mat_diag_
    );
    read_dependent_data_.numerator.mother_probability = (
      read_dependent_data_.mother_vec * somatic_probability_mat_diag_
    );
    read_dependent_data_.numerator.father_probability = (
      read_dependent_data_.father_vec * somatic_probability_mat_diag_
    );
  }
}

/**
 * Calculates denominator, probability of the observed data or numerator,
 * probability of no mutation.
 *
 * @param  is_numerator True if calculating probability of numerator.
 */
void TrioModel::GermlineTransition(bool is_numerator) {
  if (!is_numerator) {
    read_dependent_data_.denominator.child_germline_probability = (
      read_dependent_data_.denominator.child_probability *
      germline_probability_mat_
    );
    read_dependent_data_.denominator.parent_probability = KroneckerProduct(
      read_dependent_data_.denominator.mother_probability,
      read_dependent_data_.denominator.father_probability
    );
    read_dependent_data_.denominator.root_mat = TrioModel::GetRootMat(
      read_dependent_data_.denominator.child_germline_probability,
      read_dependent_data_.denominator.parent_probability
    );
    read_dependent_data_.denominator.sum = read_dependent_data_.denominator.root_mat.sum();
  } else {
    read_dependent_data_.numerator.child_germline_probability = (
      read_dependent_data_.numerator.child_probability *
      germline_probability_mat_num_
    );
    read_dependent_data_.numerator.parent_probability = KroneckerProduct(
      read_dependent_data_.numerator.mother_probability,
      read_dependent_data_.numerator.father_probability
    );
    read_dependent_data_.numerator.root_mat = TrioModel::GetRootMat(
      read_dependent_data_.numerator.child_germline_probability,
      read_dependent_data_.numerator.parent_probability
    );
    read_dependent_data_.numerator.sum = read_dependent_data_.numerator.root_mat.sum();
  }
}

/**
 * Returns root matrix for use in MutationProbability.
 *
 * @param  child_germline_probability Matrix after child germline transition.
 * @param  parent_probability         Kronecker product of both parent vectors.
 * @return                            1 x 256 final matrix at the root of tree.
 */
RowVector256d TrioModel::GetRootMat(const RowVector256d &child_germline_probability,
                                    const RowVector256d &parent_probability) {
  return child_germline_probability.cwiseProduct(
    parent_probability).cwiseProduct(population_priors_);
}

/**
 * Generates a 16 x 4 alpha frequencies matrix given the sequencing error rate
 * and dirichlet dispersion. The order of the alpha frequencies correspond to
 * the genotypes. Each alpha should sum to 1.
 * 
 * Current values are placeholders until they are estimated in Spring 2014.
 *
 * @return  16 x 4 Eigen matrix of Dirichlet multinomial alpha parameters
 *          alpha = (alpha_1, ..., alpha_K) for a K-category Dirichlet
 *          distribution (where K = 4 = kNucleotideCount) that vary with each
 *          combination of parental genotype and reference nucleotide.
 */
Matrix16_4d TrioModel::Alphas() {
  Matrix16_4d alphas;
  double homozygous = 1.0 - sequencing_error_rate_;
  double other = sequencing_error_rate_ / 3.0;
  double heterozygous = 0.5 - other;

  //        A             C             G             T
  alphas << homozygous,   other,        other,        other,
            heterozygous, heterozygous, other,        other,
            heterozygous, other,        heterozygous, other,
            heterozygous, other,        other,        heterozygous,

            heterozygous, heterozygous, other,        other,
            other,        homozygous,   other,        other,
            other,        heterozygous, heterozygous, other,
            other,        heterozygous, other,        heterozygous,
            
            heterozygous, other,        heterozygous, other,
            other,        heterozygous, heterozygous, other,
            other,        other,        homozygous,   other,
            other,        other,        heterozygous, heterozygous,

            heterozygous, other,        other,        heterozygous,
            other,        heterozygous, other,        heterozygous,
            other,        other,        heterozygous, heterozygous,
            other,        other,        other,        homozygous;
  return alphas;
}

/**
 * Returns true if the two TrioModel objects are equal to each other within
 * epsilon precision.
 * 
 * @param  other TrioModel object to be compared.
 * @return       True if the two TrioModel objects are equal to each other.
 */
bool TrioModel::Equals(const TrioModel &other) {
  bool attr_table[11] = {
    Equal(population_mutation_rate_, other.population_mutation_rate_),
    Equal(germline_mutation_rate_, other.germline_mutation_rate_),
    Equal(somatic_mutation_rate_, other.germline_mutation_rate_),
    Equal(sequencing_error_rate_, other.sequencing_error_rate_),
    Equal(dirichlet_dispersion_, other.dirichlet_dispersion_),
    nucleotide_frequencies_.isApprox(other.nucleotide_frequencies_, kEpsilon),
    genotype_mat_.isApprox(other.genotype_mat_, kEpsilon),
    population_priors_.isApprox(other.population_priors_, kEpsilon),
    germline_probability_mat_.isApprox(other.germline_probability_mat_, kEpsilon),
    somatic_probability_mat_.isApprox(other.somatic_probability_mat_, kEpsilon),
    read_dependent_data_.sequencing_probability_mat.isApprox(
      other.read_dependent_data_.sequencing_probability_mat,
      kEpsilon
    )
  };
  if (all_of(begin(attr_table), end(attr_table), [](bool i) { return i; })) {
    return true;
  } else {
    return false;
  }
}

double TrioModel::population_mutation_rate() {
  return population_mutation_rate_;
}

/**
 * Sets population_mutation_rate_, genotype_mat_ and population_priors_.
 */
void TrioModel::set_population_mutation_rate(double rate) {
  population_mutation_rate_ = rate;
  population_priors_ = PopulationPriors();
}

double TrioModel::germline_mutation_rate() {
  return germline_mutation_rate_;
}

/**
 * Sets germline_mutation_rate_, germline_probability_mat_ and
 * germline_probability_mat_num_.
 */
void TrioModel::set_germline_mutation_rate(double rate) {
  germline_mutation_rate_ = rate;
  germline_probability_mat_ = TrioModel::GermlineProbabilityMat();
  germline_probability_mat_num_ = TrioModel::GermlineProbabilityMat(true);
}

double TrioModel::somatic_mutation_rate() {
  return somatic_mutation_rate_;
}

/**
 * Sets somatic_mutation_rate_, somatic_probability_mat_ and
 * somatic_probability_mat_diag_.
 */
void TrioModel::set_somatic_mutation_rate(double rate) {
  somatic_mutation_rate_ = rate;
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMat();
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMatDiag();
}

double TrioModel::sequencing_error_rate() {
  return sequencing_error_rate_;
}

/**
 * Sets sequencing_error_rate_ and alphas_.
 */
void TrioModel::set_sequencing_error_rate(double rate) {
  sequencing_error_rate_ = rate;
  alphas_ = TrioModel::Alphas();
}

double TrioModel::dirichlet_dispersion() {
  return dirichlet_dispersion_;
}

/**
 * Sets dirichlet_dispersion_ and alphas_.
 */
void TrioModel::set_dirichlet_dispersion(double dispersion) {
  dirichlet_dispersion_ = dispersion;
  alphas_ = TrioModel::Alphas();
}

RowVector4d TrioModel::nucleotide_frequencies() {
  return nucleotide_frequencies_;
}

/**
 * Sets nucleotide_frequencies_ and population_priors_.
 */
void TrioModel::set_nucleotide_frequencies(const RowVector4d &frequencies) {
  nucleotide_frequencies_ = frequencies;
  population_priors_ = TrioModel::PopulationPriors();
}

bool TrioModel::has_mutation() {
  return read_dependent_data_.has_mutation;
}

void TrioModel::set_has_mutation(bool has_mutation) {
  read_dependent_data_.has_mutation = has_mutation;
}

RowVector16d TrioModel::genotype_mat() {
  return genotype_mat_;
}

RowVector256d TrioModel::population_priors() {
  return population_priors_;
}

Matrix16_256d TrioModel::germline_probability_mat() {
  return germline_probability_mat_;
}

Matrix16_16d TrioModel::somatic_probability_mat() {
  return somatic_probability_mat_;
}

Matrix3_16d TrioModel::sequencing_probability_mat() {
  return read_dependent_data_.sequencing_probability_mat;
}

Matrix16_4d TrioModel::alphas() {
  return alphas_;
}
