/**
 * @file simulation_driver.cc
 * @author Melissa Ip
 *
 * This file runs simulations using command line inputs and the SimulationModel
 * class. The output file is formatted into two columns. Probability of mutation
 * is the first column and whether the site has a mutation is the second column.
 * Each entry is tab deliminated and each site is placed on a new line. For
 * example:
 *
 * 1.6732e-10      0
 * 0.00709331      1
 *
 * To compile on Herschel, use the following command to include the GSL library:
 * c++ -std=c++11  -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o simulation_driver simulation_driver.cc
 *
 * To run this file, provide the following command line inputs:
 * ./simulation_driver <output.txt> <#samples> <coverage> <germline mutation rate> <somatic mutation rate>
 */
#include <cerrno>
#include <stdio.h>
#include <time.h>

#include "simulation_model.cc"


/**
 * Prints user defined errors and exits the program.
 *
 * @param  msg Message to be printed.
 */
void Die(const char *msg) {
  if (errno == EDOM) {
    perror(msg);
  } else {
    printf("ERROR: %s\n", msg);
  }
  exit(EXIT_FAILURE);
}

int main(int argc, const char *argv[]) {
  if (argc < 6) {        
    Die("USAGE: simulation_driver <output.txt> <#samples> <coverage> "
        "<germline mutation rate> <somatic mutation rate>");
  }

  const string file_name = argv[1];
  const int experiment_count = atoi(argv[2]);
  const unsigned int coverage = strtoull(argv[3], NULL, 10);
  const double germline_mutation_rate = strtod(argv[4], NULL);
  const double somatic_mutation_rate = strtod(argv[5], NULL);

  srand(time(NULL));
  SimulationModel sim(coverage, germline_mutation_rate, somatic_mutation_rate);
  sim.WriteProbability(file_name, experiment_count);

  return 0;
}
