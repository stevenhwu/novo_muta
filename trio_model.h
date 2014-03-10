/**
 * @file trio_model.h
 * @author Melissa Ip
 *
 * The TrioModel class is used to create an object that contains sequencing
 * read data, estimation parameters, and useful matrices in calculating the
 * probability of de novo mutation.
 *
 * The probability of mutation is calculated using a modified trio model,
 * similar but different to the model described in the following paper:
 * 
 * Cartwright et al.: Family-Based Method for Capturing De Novo Mutations
 * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3728889/
 *
 * Example usage:
 *
 *   TrioModel params;
 *   ReadDataVector data = {
 *     {30, 0, 0, 0},
 *     {30, 0, 0, 0},
 *     {30, 0, 0, 0}
 *   };
 *   double probability = params.MutationProbability(data);
 *   data.set_germline_mutation_rate(0.000001);
 *   double new_probability = params.MutationProbability(data);
 */
#ifndef TRIO_MODEL_H
#define TRIO_MODEL_H

 
#include "utilities.cc"

/**
 * TrioModel class header. See top of file for a complete description.
 */
class TrioModel {
public:
  // Default constructor and constructor to customize parameters.
  TrioModel();
  TrioModel(double population_mutation_rate, double germline_mutation_rate,
            double somatic_mutation_rate, double sequencing_error_rate,
            double dirichlet_dispersion,
            const RowVector4d &nucleotide_frequencies);

  // Calculates probability of mutation given input data.
  double MutationProbability(const ReadDataVector &data_vec);

  // True if the two TrioModel objects are equal to each other.
  bool Equals(const TrioModel &other);

    // Get and set methods.
  double population_mutation_rate();
  void set_population_mutation_rate(double rate);
  double germline_mutation_rate();
  void set_germline_mutation_rate(double rate);
  double somatic_mutation_rate();
  void set_somatic_mutation_rate(double rate);
  double sequencing_error_rate();
  void set_sequencing_error_rate(double rate);
  double dirichlet_dispersion();
  void set_dirichlet_dispersion(double dispersion);
  RowVector4d nucleotide_frequencies();
  void set_nucleotide_frequencies(const RowVector4d &frequencies);
  bool has_mutation();
  void set_has_mutation(bool has_mutation);
  RowVector16d genotype_mat();
  RowVector256d population_priors();
  Matrix16_256d germline_probability_mat();
  Matrix16_16d somatic_probability_mat();
  Matrix3_16d sequencing_probability_mat();
  Matrix16_4d alphas();

  // E-step methods.
  RowVector16d GetSomaticStatistic();

private:
  // Helper methods for MutationProbability.
  void GermlineTransition(bool is_numerator=false);
  void SomaticTransition(bool is_numerator=false);
  RowVector256d GetRootMat(const RowVector256d &child_germline_probability,
                           const RowVector256d &parent_probability);

  // Calculates probability of allele spectrum given read counts.
  double SpectrumProbability(const RowVector4d &nucleotide_counts);

  // Methods for setting up the model and relevant arrays.
  RowVector256d PopulationPriors();
  double GermlineMutation(int child_nucleotide_idx, int parent_genotype_idx,
                          bool no_mutation_flag);
  Matrix16_256d GermlineProbabilityMat(bool no_mutation_flag=false);
  double SomaticMutation(int nucleotide_idx, int other_nucleotide_idx);
  Matrix16_16d SomaticProbabilityMat();
  Matrix16_16d SomaticProbabilityMatDiag();
  void SequencingProbabilityMat(const ReadDataVector &data_vec);
  Matrix16_4d Alphas();
  
  // E-step methods.
  double SomaticMutationCountsMatrix();

  // Instance variables.
  double population_mutation_rate_;
  double germline_mutation_rate_;
  double somatic_mutation_rate_;
  double sequencing_error_rate_;
  double dirichlet_dispersion_;  // unused
  RowVector4d nucleotide_frequencies_;
  Matrix16_4d alphas_;
  RowVector16d genotype_mat_;  // unused
  RowVector256d population_priors_;
  Matrix16_256d germline_probability_mat_;
  Matrix16_256d germline_probability_mat_num_;
  Matrix16_16d somatic_probability_mat_;
  Matrix16_16d somatic_probability_mat_diag_;

  // TODO: Update simulation files to use refactored trio model.
  struct ReadDependentData {
    vector<double> max_elements;  // stores max element of sequencing_probability_mat when rescaling to normal space
    Matrix3_16d sequencing_probability_mat;  // P(R|somatic genotype)
    RowVector16d child_vec;
    RowVector16d mother_vec;
    RowVector16d father_vec;
    struct TreePeels {
      RowVector16d child_probability;  // P(R|zygotic genotype)
      RowVector16d mother_probability;
      RowVector16d father_probability;
      RowVector256d child_germline_probability;
      RowVector256d parent_probability;  // P(R|mom and dad genotype)
      RowVector256d root_mat;  // P(R|mom and dad genotype) * P(mom and dad genotype)
      double sum;  // P(R)
    } denominator, numerator;
    bool has_mutation;  // simulation only
  } read_dependent_data_;
};

#endif
