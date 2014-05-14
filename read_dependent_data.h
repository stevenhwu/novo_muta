/**
 * @file read_dependent_data.h
 * @author Melissa Ip
 *
 * The ReadDependentData class contains the members that are dependent on
 * ReadData/ReadDataVector including matrices and vectors that are calculated
 * from sequencing error and by the tree-peeling algorithm. 
 *
 * TODO: The simulation program needs to be refactored to use this class.
 */
#ifndef READ_DEPENDENT_DATA_H
#define READ_DEPENDENT_DATA_H

#include "utility.h"


/**
 * ReadDependentData class header. See top of file for a complete description.
 */
class ReadDependentData {
 public:
  ReadDependentData();  // Default constructor leaves read_data_vec empty.
  ReadDependentData(const ReadDataVector &data_vec);  // Constructor that initializes read_data_vec.

  // Instance member variables.
  ReadDataVector read_data_vec;
  vector<double> max_elements;  // Stores max element of sequencing_probability_mat when rescaling to normal space.
  Matrix3_16d sequencing_probability_mat;  // P(R|somatic genotype)
  RowVector16d child_somatic_probability;
  RowVector16d mother_somatic_probability;
  RowVector16d father_somatic_probability;
  struct TreePeel {  // TODO: Change to TreePeel class.
    RowVector16d child_zygotic_probability; // P(R|zygotic genotype)
    RowVector16d mother_zygotic_probability;
    RowVector16d father_zygotic_probability;
    RowVector256d child_germline_probability;
    RowVector256d parent_probability; // P(R|mom and dad genotype)
    RowVector256d root_mat; // P(R|mom and dad genotype) * P(mom and dad genotype)
    double sum; // P(R)
  } denominator, numerator;
  bool has_mutation;  // Simulation only.
};

#endif
