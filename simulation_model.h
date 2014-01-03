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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <fstream>

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
  
  // Generates random samples and probabilities in text file.
  void WriteProbability(const string &file_name, int experiment_count);

private:
  // Simulation methods.
  int Mutate(int genotype_idx, bool is_germline=false,
             int parent_genotype_idx=-1);
  ReadData DirichletMultinomialSample(int genotype_idx);

  // Instance variables.
  TrioModel params_;  // Default initialization.
  unsigned int coverage_;
  bool has_mutation_;
};
