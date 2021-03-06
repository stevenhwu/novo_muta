/**
 * @file pileup_driver.cc
 * @author Melissa Ip
 *
 * This file accepts 3 input pileup files. The pileup data is read and parsed
 * into sequencing read data that the TrioModel can process. Default parameter
 * values are used. Assume that all pileup files have the same number of sites
 * and thus can be aligned.
 * 
 * To compile on Herschel without using cmake and include GSL:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -I/usr/local/include -o pileup_driver pileup_utility.cc pileup_driver.cc
 *
 * To run this file, provide the following command line inputs:
 * ./pileup_driver <output>.txt <child>.pileup <mother>.pileup <father>.pileup
 *
 * See top of pileup_utility.h for additional information.
 */
#include "pileup_utility.h"


int main(int argc, const char *argv[]) {
  if (argc < 5) {
    Die("USAGE: pileup_driver <output>.txt <child>.pileup <mother>.pileup "
        "<father>.pileup");
  }

  const string file_name = argv[1];
  const string child_pileup = argv[2];
  const string mother_pileup = argv[3];
  const string father_pileup = argv[4];

  ProcessPileup(file_name, child_pileup, mother_pileup, father_pileup);

  return 0;
}
