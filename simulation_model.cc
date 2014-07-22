/**
 * @file simulation_model.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the SimulationModel class.
 *
 * See top of simulation_model.h for a complete description.
 */
#include "simulation_model.h"


/**
 * Constructor to initialize SimulationModel object that contains necessary
 * parameters for the random generation of samples and probabilities. Currently
 * only coverage and germline and somatic mutation rates can be modified. All
 * other parameters are set to default.
 *
 * has_mutation_ keeps track of whether each site contains a mutation and is
 * reused for all sites.
 *
 * @param  coverage               Coverage.
 * @param  germline_mutation_rate Germline mutation rate.
 * @param  somatic_mutation_rate  Somatic mutation rate.
 */
SimulationModel::SimulationModel(unsigned int coverage,
                                 double germline_mutation_rate,
                                 double somatic_mutation_rate)
    :  coverage_{coverage}, has_mutation_{false} {
  params_.set_germline_mutation_rate(germline_mutation_rate);
  params_.set_somatic_mutation_rate(somatic_mutation_rate);
}

/**
 * Seeds rand function with time and GSL functions with random number generator.
 */
void SimulationModel::Seed() {
  srand(time(NULL));
  // Creates a generator chosen by the environment variable GSL_RNG_TYPE.
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

/**
 * Frees GSL random number generator.
 */
void SimulationModel::Free() {
  gsl_rng_free(r);
}

/**
 * Mutates a numeric genotype based on either germline or somatic transition
 * matrix in the TrioModel.
 *
 * is_germline is false by default to process somatic mutation. If this is set
 * to true, then this method will process germline mutation and assume
 * parent_genotype_idx is not -1.
 *
 * @param  genotype_idx        Index of genotype.
 * @param  is_germline         False by default. Set to true to process germline
 *                             mutation.
 * @param  parent_genotype_idx Index of parent genotype.
 * @return                     Index of mutated genotype.
 */
int SimulationModel::Mutate(int genotype_idx, bool is_germline,
                            int parent_genotype_idx) {
  // Sets probability matrices to use either germline or somatic probabilities.
  RowVector16d mat = RowVector16d::Zero();
  if (!is_germline) {
    mat = params_.somatic_probability_mat().row(genotype_idx);
  } else {
    mat = params_.germline_probability_mat().col(parent_genotype_idx);
  }

  // Randomly mutates the genotype using the probabilities as weights.
  int mutated_genotype_idx = SimulationModel::RandomDiscreteChoice(
    kGenotypeCount,
    mat
  );
  if (mutated_genotype_idx != genotype_idx) {
    has_mutation_ = true;
  }
  return mutated_genotype_idx;
}

/**
 * Generates a numeric child genotype by picking a random allele from each of
 * the given numeric parent genotypes.
 *
 * @param  mother_genotype Numeric mother genotype.
 * @param  father_genotype Numeric father genotype.
 * @return                 Numeric child genotype.
 */
int SimulationModel::GetChildGenotype(int mother_genotype, int father_genotype) {
  int child_allele1 = kGenotypeNumIndex(mother_genotype, rand() % 2);
  int child_allele2 = kGenotypeNumIndex(father_genotype, rand() % 2);
  for (int i = 0; i < kGenotypeCount; ++i) {
    if (child_allele1 == i / kNucleotideCount &&
        child_allele2 == i % kNucleotideCount) {
      return i;
    }
  }
  return -1;  // ERROR: This should not happen.
}

/**
 * Uses alpha frequencies based on the somatic genotype to select nucleotide
 * frequencies and uses these frequencies to draw sequencing reads at a
 * specified coverage (Dirichlet multinomial). K is kNucleotideCount.
 *
 * @param  genotype_idx Index of genotype.
 * @return              Read counts drawn from Dirichlet multinomial.
 */
ReadData SimulationModel::DirichletMultinomialSample(int genotype_idx) {
  // Converts alpha to double array.
  auto alpha_vec = params_.alphas().row(genotype_idx);
  const double alpha[kNucleotideCount] = {alpha_vec(0), alpha_vec(1),
                                          alpha_vec(2), alpha_vec(3)};

  // Sets alpha frequencies using dirichlet distribution in theta.
  double theta[kNucleotideCount] = {0.0};
  gsl_ran_dirichlet(r, kNucleotideCount, alpha, theta);

  // Sets sequencing reads using multinomial distribution in reads.
  unsigned int reads[kNucleotideCount] = {0};
  gsl_ran_multinomial(r, kNucleotideCount, coverage_, theta, reads);

  // Converts reads to ReadData.
  ReadData data = {0};
  for (int i = 0; i < kNucleotideCount; ++i) {
    data.reads[i] = reads[i];
  }
  return data;
}

/**
 * Generates a 3 x size Eigen matrix of random genotypes for child, mother,
 * and father.
 *
 * @param  size Number of genotypes to generate.
 * @return      3 x size Eigen matrix of genotypes.
 */
MatrixXi SimulationModel::GetGenotypesMatrix(int size) {
  // Generates random samples using population priors as weights.
  RowVectorXi parent_genotypes = SimulationModel::RandomDiscreteChoice(
    kGenotypePairCount,
    params_.population_priors(),
    size
  );

  // Extracts parent genotypes from samples and gets child genotypes.
  MatrixXi genotypes_mat(3, size);
  genotypes_mat.row(1) = parent_genotypes / kGenotypeCount;
  for (int i = 0; i < size; ++i) {
    int father_genotype = parent_genotypes(i) % kGenotypeCount;
    genotypes_mat.row(2)(i) = father_genotype;
    genotypes_mat.row(0)(i) = SimulationModel::GetChildGenotype(
      genotypes_mat(1, i),  // Mother genotype.
      father_genotype
    );
  }

  return genotypes_mat;
}

/**
 * Generates size random trios and keeps track of whether each ReadDataVector
 * has a mutation by adding it to the has_mutation_vec_.
 *
 * @param  size Number of random trios.
 * @return      TrioVector containing random trios.
 */
TrioVector SimulationModel::GetRandomTrios(int size) {
  TrioVector random_trios;
  MatrixXi genotypes_mat = SimulationModel::GetGenotypesMatrix(size);
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);

  for (int i = 0; i < size; ++i) {
    int child_genotype = genotypes_mat(0, i);
    int mother_genotype = genotypes_mat(1, i);
    int father_genotype = genotypes_mat(2, i);

    // Processes germline mutation. Germline matrix requires no Kronecker.
    int child_germline_genotype = SimulationModel::Mutate(
      child_genotype,
      true,
      mother_genotype*kGenotypeCount + father_genotype
    );

    // Processes somatic mutation.
    int child_somatic_genotype = SimulationModel::Mutate(child_germline_genotype);
    int mother_somatic_genotype = SimulationModel::Mutate(mother_genotype);
    int father_somatic_genotype = SimulationModel::Mutate(father_genotype);

    // Creates reads from somatic genotypes using the Dirichlet multinomial.
    ReadData child_read = SimulationModel::DirichletMultinomialSample(child_somatic_genotype);
    ReadData mother_read = SimulationModel::DirichletMultinomialSample(mother_somatic_genotype);
    ReadData father_read = SimulationModel::DirichletMultinomialSample(father_somatic_genotype);
    ReadDataVector data_vec = {child_read, mother_read, father_read};
    random_trios.push_back(data_vec);

    // Records has_mutation_ in order relevant vector.
    int trio_index = IndexOfReadDataVector(data_vec, trio_vec);
    if (trio_index != -1) {
      mutation_table_[trio_index].push_back(has_mutation_);
      has_mutation_vec_.push_back(has_mutation_);  // Must be valid trio for record.
    }
    has_mutation_ = false;  // Resets for the next simulation.
  }
  return random_trios;
}

/**
 * Generates size random samples using population priors as weights and outputs
 * their probabilities and whether that site contains a mutation
 * (1=true, 0=false) to a text file. The file is tab separated and each site
 * is on a new line.
 *
 * @param  file_name File name.
 * @param  size      Number of experiments or trios.
 */
void SimulationModel::WriteProbability(const string &file_name, int size) {
  ofstream fout(file_name);
  TrioVector random_trios = SimulationModel::GetRandomTrios(size);
  for (int i = 0; i < size; ++i) {
    double probability = params_.MutationProbability(random_trios[i]);
    fout << probability << "\t" << has_mutation_vec_[i] << "\n";
  }
  fout.close();
}

/**
 * Writes to a text file the index of the key trio, how many random trios had a
 * mutation, how many random trios had no mutation, tab separated, each trio
 * is placed on a new line. Assumes 4x coverage.
 *
 * @param  file_name File name.
 * @param  size      Number of random trios.
 */
void SimulationModel::WriteMutationCounts(const string &file_name, int size) {
  ofstream fout(file_name);
  TrioVector random_trios = SimulationModel::GetRandomTrios(size);
  for (int i = 0; i < kTrioCount; ++i) {
    vector<bool> mutations = mutation_table_[i];
    int has_mutation_total = 0;
    int has_no_mutation_total = 0;
    if (mutations.size() > 0) {
      has_mutation_total = count(mutations.begin(), mutations.end(), true);
      has_no_mutation_total = mutations.size() - has_mutation_total;
    }

    fout << i << "\t" << has_mutation_total << "\t"
         << has_no_mutation_total << "\n";
  }
  fout.close();
}

/**
 * Writes to stdout the index of the key trio, how many random trios had a
 * mutation, how many random trios had no mutation, tab separated, each trio
 * is placed on a new line. Assumes 4x coverage.
 *
 * @param  size Number of random trios.
 */
void SimulationModel::PrintMutationCounts(int size) {
  TrioVector random_trios = SimulationModel::GetRandomTrios(size);
  for (int i = 0; i < kTrioCount; ++i) {
    vector<bool> mutations = mutation_table_[i];
    int has_mutation_total = 0;
    int has_no_mutation_total = 0;
    if (mutations.size() > 0) {
      has_mutation_total = count(mutations.begin(), mutations.end(), true);
      has_no_mutation_total = mutations.size() - has_mutation_total;
    }

    cout << i << "\t" << has_mutation_total << "\t"
         << has_no_mutation_total << "\n";
  }
}

/**
 * Generates a random sample from the range [0, K) using K size probabilities
 * as weights.
 *
 * @param  K              Number of cateogories in probabilities.
 * @param  probabilitites RowVector of probabilities associated with each entry
 *                        in the samples generated by K.
 * @return                Random element based on probabilities from samples
 *                        generated by K.
 */
int SimulationModel::RandomDiscreteChoice(size_t K,
                                          const RowVectorXd &probabilities) {
  // Converts probabilities to double array p.
  int length = probabilities.size();
  double p[length];
  for (int i = 0; i < length; ++i) {
    p[i] = probabilities(i);
  }

  // Creates discrete random number generator.
  gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(K, p);
  size_t random = gsl_ran_discrete(r, g);
  gsl_ran_discrete_free(g);
  
  return (int) random;
}

/**
 * Generates RowVector of random samples from the range [0, K) using K size
 * probabilities as weights.
 *
 * @param  K              Number of cateogories in probabilities.
 * @param  probabilitites RowVector of probabilities associated with each entry
 *                        in the samples generated by K.
 * @param  size           Number of samples.
 * @return                Random samples based on probabilities from samples
 *                        generated by K.
 */
RowVectorXi SimulationModel::RandomDiscreteChoice(
    size_t K, const RowVectorXd &probabilities, int size) {
  // Creates size RowVector to hold random samples.
  RowVectorXi random_samples(size);
  for (int i = 0; i < size; ++i) {
    random_samples(i) = SimulationModel::RandomDiscreteChoice(K, probabilities);
  }
  return random_samples;
}

unsigned int SimulationModel::coverage() {
  return coverage_;
}

void SimulationModel::set_coverage(unsigned int coverage) {
  coverage_ = coverage;
}

double SimulationModel::germline_mutation_rate() {
  return params_.germline_mutation_rate();
}

void SimulationModel::set_germline_mutation_rate(double rate) {
  params_.set_germline_mutation_rate(rate);
}

double SimulationModel::somatic_mutation_rate() {
  return params_.somatic_mutation_rate();
}

void SimulationModel::set_somatic_mutation_rate(double rate) {
  params_.set_somatic_mutation_rate(rate);
}

bool SimulationModel::has_mutation() {
  return has_mutation_;
}

void SimulationModel::set_has_mutation(bool has_mutation) {
  has_mutation_ = has_mutation;
}
