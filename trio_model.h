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
            double dirichlet_dispersion, const Array4d &nucleotide_frequencies);

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
  Array4d nucleotide_frequencies();
  void set_nucleotide_frequencies(const Array4d &frequencies);
  Array16d genotype_mat();
  Array256d population_priors();
  Array16_256d germline_probability_mat();
  Array16_16d somatic_probability_mat();
  Array3_16d sequencing_probability_mat();

private:
  // Methods for setting up the model and relevant arrays.
  Array256d PopulationPriors();
  double GermlineMutation(int child_nucleotide_idx, int parent_genotype_idx,
                          bool no_mutation_flag);
  Array16_256d GermlineProbabilityMat(bool no_mutation_flag=false);
  double SomaticMutation(int nucleotide_idx, int other_nucleotide_idx);
  Array16_16d SomaticProbabilityMat();
  Array3_16d SequencingProbabilityMat(const ReadDataVector &data_vec);

  // Instance variables.
  double population_mutation_rate_;
  double germline_mutation_rate_;
  double somatic_mutation_rate_;
  double sequencing_error_rate_;
  double dirichlet_dispersion_;
  Array4d nucleotide_frequencies_;
  vector<double> max_elements_;
  Array16d genotype_mat_;
  Array256d population_priors_;
  Array16_256d germline_probability_mat_;
  Array16_16d somatic_probability_mat_;
  Array3_16d sequencing_probability_mat_;
};
