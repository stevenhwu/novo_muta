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
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMat();
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
                     const Array4d &nucleotide_frequencies)
    : population_mutation_rate_{population_mutation_rate},
      germline_mutation_rate_{germline_mutation_rate},
      somatic_mutation_rate_{somatic_mutation_rate},
      sequencing_error_rate_{sequencing_error_rate},
      dirichlet_dispersion_{dirichlet_dispersion},
      nucleotide_frequencies_{nucleotide_frequencies} {
  population_priors_ = TrioModel::PopulationPriors();
  germline_probability_mat_ = TrioModel::GermlineProbabilityMat();
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMat();
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
 * Initializes and updates sequencing_probability_mat_ using
 * sequencing_error_rate_ and dirichlet_dispersion_.
 *
 * @param   data_vec Read counts in order of child, mother and father.
 * @return           Probability of mutation given read data and parameters.
 */
double TrioModel::MutationProbability(const ReadDataVector &data_vec) {
  // Creates and sets sequencing probabilities given read data.
  sequencing_probability_mat_ = TrioModel::SequencingProbabilityMat(data_vec);
  Array16d child_vec = sequencing_probability_mat_.row(0);
  Array16d mother_vec = sequencing_probability_mat_.row(1);
  Array16d father_vec = sequencing_probability_mat_.row(2);

  // Multiplies vectors by somatic transition matrix.
  Array16d child_probability = DotProduct(child_vec, somatic_probability_mat_);
  Array16d mother_probability = DotProduct(mother_vec, somatic_probability_mat_);
  Array16d father_probability = DotProduct(father_vec, somatic_probability_mat_);

  // Calculates denominator, probability of the observed data.
  Array256d child_germline_probability = DotProduct(
    child_probability,
    germline_probability_mat_
  );
  Array256d parent_probability = KroneckerProduct(
    mother_probability,
    father_probability
  );
  Array256d demoninator_mat = (child_germline_probability * parent_probability *
    population_priors_);
  double demoninator_sum = demoninator_mat.sum();

  // Calculates numerator, probability of no mutation.
  Array16_16d somatic_probability_diag = GetDiagonal(somatic_probability_mat_);
  Array16d child_probability_num = DotProduct(child_vec,
                                              somatic_probability_diag);
  Array16d mother_probability_num = DotProduct(mother_vec,
                                               somatic_probability_diag);
  Array16d father_probability_num = DotProduct(father_vec,
                                               somatic_probability_diag);

  Array16_256d germline_probability_mat_num = TrioModel::GermlineProbabilityMat(true);
  Array256d child_germline_probability_num = DotProduct(
    child_probability_num,
    germline_probability_mat_num
  );
  Array256d parent_probability_num = KroneckerProduct(
    mother_probability_num,
    father_probability_num
  );
  Array256d numerator_mat = (child_germline_probability_num *
    parent_probability_num * population_priors_);
  double numerator_sum = numerator_mat.sum();
  
  return 1 - numerator_sum/demoninator_sum;
}

/**
 * Returns 1 x 256 Eigen probability array. This is an order-relevant
 * representation of the possible events in the sample space that covers all
 * possible parent genotype combinations. For example:
 *
 * [P(AAAA), P(AAAC), P(AAAG), P(AAAT), P(AACA), P(AACC), P(AACG)...]
 *
 * Sets genotype_mat_ as the 1 x 16 Eigen probability array for a single parent.
 *
 * Creates nucleotide mutation frequencies {alpha_A, alpha_C, alpha_G, alpha_T}
 * based on the nucleotide frequencies and population mutation rate (theta).
 * These frequencies and nucleotide counts {n_A, n_C, n_G, n_T} are used in the
 * Dirichlet multinomial.
 *
 * For example, both parents have genotype AA, resulting in N = 4:
 *
 * population_mutation_rate_ = 0.00025;
 * nucleotide_frequencies_ << 0.25, 0.25, 0.25, 0.25;
 * nucleotide_counts = {4, 0, 0, 0};
 *
 * @return  1 x 256 Eigen probability array in log e space where the (i, j)
 *          element is the probability that the mother has genotype i and the
 *          father has genotype j.
 */
Array256d TrioModel::PopulationPriors() {
  // Calculates nucleotide mutation frequencies using given mutation rate.
  Array4d nucleotide_mutation_frequencies = (nucleotide_frequencies_ *
    population_mutation_rate_);

  Array16_16d population_priors = Array16_16d::Zero();
  // Resizes population_priors to 1 x 256 array.
  Array256d population_priors_flattened;
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      // Convert nucleotide_counts to ReadData for DirichletMultinomialLog().
      Array4d nucleotide_counts = kTwoParentCounts(i, j);
      ReadData nucleotide_read;
      for (int k = 0; k < kNucleotideCount; ++k) {
        nucleotide_read.reads[k] = nucleotide_counts(k);
      }
      // Calculates probability using the Dirichlet multinomial in normal space.
      double log_probability = DirichletMultinomialLog(
        nucleotide_mutation_frequencies,
        nucleotide_read
      );
      population_priors(i, j) = exp(log_probability);
      int idx = i * kGenotypeCount + j;
      population_priors_flattened(idx) = population_priors(i, j);
    }
  }

  // Sets genotype_mat_ to collapsed single parent probability matrix.
  genotype_mat_ = population_priors.rowwise().sum();

  return population_priors_flattened;
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
  bool is_in_array_flag = false;
  if (child_nucleotide_idx == parent_genotype(0) ||
      child_nucleotide_idx == parent_genotype(1)) {
    is_in_array_flag = true;
  }

  // Determines if the comparison is homozygous, heterozygous or no match.
  if (is_in_array_flag) {
    if (parent_genotype(0) == parent_genotype(1)) {
      return homozygous_match;
    } else {
      if (no_mutation_flag) {
        return homozygous_match/2;
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
 * Calculates the probability array for the offspring using the given mutation
 * rate, derived from the Kronecker product of a 4 x 16 Eigen parent 
 * probability array with itself.
 *
 * This is a transition matrix used to mutate the child germline genotype.
 *
 * @param  no_mutation_flag False by default. Set to true to calculate
 *                          probability of no mutation.
 * @return                  16 x 256 Eigen probability array.
 */
Array16_256d TrioModel::GermlineProbabilityMat(bool no_mutation_flag) {
  Array4_16d germline_probability_mat = Array4_16d::Zero();
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
 * of each other and calculates somatic mutation probability array for a 
 * single parent.
 *
 * This is a transition matrix used to mutate the somatic genotypes.
 *
 * @return  16 x 16 Eigen probability array where the first dimension is the
 *          original somatic genotypes and the second dimension is the mutated 
 *          genotype.
 */
Array16_16d TrioModel::SomaticProbabilityMat() {
  Array44d somatic_probability_mat = Array44d::Zero();
  for (int i = 0; i < kNucleotideCount; ++i) {
    for (int j = 0; j < kNucleotideCount; ++j) {
      double probability = TrioModel::SomaticMutation(i, j);
      somatic_probability_mat(i, j) = probability;
    }
  }
  return KroneckerProduct(somatic_probability_mat);
}

/**
 * Calculates the probability of sequencing error for all read data. Assume
 * data contains 3 reads (child, mother, father). Assume each chromosome is
 * equally likely to be sequenced.
 *
 * Adds the max element of all reads in ReadDataVector to max_elements_ before
 * rescaling to normal space.
 *
 * @param  data_vec ReadDataVector containing nucleotide counts for trio family.
 * @return          3 x 16 Eigen probability array.
 */
Array3_16d TrioModel::SequencingProbabilityMat(const ReadDataVector &data_vec) {
  Array3_16d sequencing_probability_mat = Array3_16d::Zero();
  Array16_4d alphas = GetAlphas(sequencing_error_rate_) * dirichlet_dispersion_;
  for (int read = 0; read < 3; ++read) {
    for (int genotype_idx = 0; genotype_idx < kGenotypeCount; ++genotype_idx) {
      auto alpha = alphas.row(genotype_idx);
      double log_probability = DirichletMultinomialLog(alpha, data_vec[read]);
      sequencing_probability_mat(read, genotype_idx) = log_probability;
    }
  }
  // Rescales to normal space and records max element of all 3 reads together.
  double max_element = sequencing_probability_mat.maxCoeff();
  max_elements_.push_back(max_element);
  return exp(sequencing_probability_mat - max_element);
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
    sequencing_probability_mat_.isApprox(other.sequencing_probability_mat_, kEpsilon)
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
 * Sets germline_mutation_rate_ and germline_probability_mat_.
 */
void TrioModel::set_germline_mutation_rate(double rate) {
  germline_mutation_rate_ = rate;
  germline_probability_mat_ = TrioModel::GermlineProbabilityMat();
}

double TrioModel::somatic_mutation_rate() {
  return somatic_mutation_rate_;
}

/**
 * Sets somatic_mutation_rate_ and somatic_probability_mat_.
 */
void TrioModel::set_somatic_mutation_rate(double rate) {
  somatic_mutation_rate_ = rate;
  somatic_probability_mat_ = TrioModel::SomaticProbabilityMat();
}

double TrioModel::sequencing_error_rate() {
  return sequencing_error_rate_;
}

void TrioModel::set_sequencing_error_rate(double rate) {
  sequencing_error_rate_ = rate;
}

Array4d TrioModel::nucleotide_frequencies() {
  return nucleotide_frequencies_;
}

/**
 * Sets nucleotide_frequencies_ and population_priors_.
 */
void TrioModel::set_nucleotide_frequencies(const Array4d &frequencies) {
  nucleotide_frequencies_ = frequencies;
  population_priors_ = TrioModel::PopulationPriors();
}

double TrioModel::dirichlet_dispersion() {
  return dirichlet_dispersion_;
}

void TrioModel::set_dirichlet_dispersion(double dispersion) {
  dirichlet_dispersion_ = dispersion;
}

Array16d TrioModel::genotype_mat() {
  return genotype_mat_;
}

Array256d TrioModel::population_priors() {
  return population_priors_;
}

Array16_256d TrioModel::germline_probability_mat() {
  return germline_probability_mat_;
}

Array16_16d TrioModel::somatic_probability_mat() {
  return somatic_probability_mat_;
}

Array3_16d TrioModel::sequencing_probability_mat() {
  return sequencing_probability_mat_;
}
