/**
 * @file simulation_model.h
 * @author Melissa Ip
 *
 * The SimulationModel class is used to validate TrioModel results using
 * default parameter values described in the TrioModel class. It generates a
 * random family pedigree based on population priors and calculates the
 * probability of mutation using the generated sample (sequencing reads are
 * drawn from the Dirichlet multinomial).
 */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <time.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "trio_model.cc"


/**
 * SimulationModel class header. See top of file for a complete description.
 */
class SimulationModel {
public:
  // No default constructor, because all parameters are given through command
  // line inputs.
  SimulationModel(unsigned int coverage, double germline_mutation_rate,
                  double somatic_mutation_rate);
  
  // Seeds random number generator during initialization.
  void Seed();
  void Free();
  // Generates random samples and probabilities in text file.
  void WriteProbability(const string &file_name, int size);

  void PrintMutationCounts(int size);
  // Get and set methods.
  unsigned int coverage();
  void set_coverage(unsigned int coverage);
  double germline_mutation_rate();
  void set_germline_mutation_rate(double rate);
  double somatic_mutation_rate();
  void set_somatic_mutation_rate(double rate);

private:
  // Simulation methods.
  int Mutate(int genotype_idx, bool is_germline=false,
             int parent_genotype_idx=-1);
  int GetChildGenotype(int mother_genotype, int father_genotype);
  ReadData DirichletMultinomialSample(int genotype_idx);
  vector<RowVectorXi> GetGenotypesMatrix(int size);
  TrioVector GetRandomTrios(int size);
  int RandomDiscreteChoice(size_t K, const RowVectorXd &probabilities);
  RowVectorXi RandomDiscreteChoice(size_t K, const RowVectorXd &probabilities,
                                   int size);

  // Instance variables.
  TrioModel params_;  // Default initialization.
  unsigned int coverage_;
  vector<bool> has_mutation_vec_;
  vector<bool> mutation_table_[42875];
  gsl_rng *r;
};
