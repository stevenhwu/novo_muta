/**
 * @file read_dependent_data.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the ReadDependentData class.
 *
 * See top of read_dependent_data.h for a complete description.
 */
#include "read_dependent_data.h"


/**
 * Default constructor.
 */
ReadDependentData::ReadDependentData()
    : has_mutation{false} {
  sequencing_probability_mat = Matrix3_16d::Zero();
  child_somatic_probability = RowVector16d::Zero();
  mother_somatic_probability = RowVector16d::Zero();
  father_somatic_probability = RowVector16d::Zero();
};

/**
 * Constructor takes in ReadDataVector.
 */
ReadDependentData::ReadDependentData(const ReadDataVector &data_vec)
    : has_mutation{false} {
  read_data_vec = data_vec;
  sequencing_probability_mat = Matrix3_16d::Zero();
  child_somatic_probability = RowVector16d::Zero();
  mother_somatic_probability = RowVector16d::Zero();
  father_somatic_probability = RowVector16d::Zero();
};
