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
ReadDependentData::ReadDependentData() {
  sequencing_probability_mat = Matrix3_16d::Zero();
  child_somatic_probability = RowVector16d::Zero();
  mother_somatic_probability = RowVector16d::Zero();
  father_somatic_probability = RowVector16d::Zero();
}

/**
 * Constructor takes in ReadDataVector.
 */
ReadDependentData::ReadDependentData(const ReadDataVector &data_vec)
    : read_data_vec{data_vec} {
  sequencing_probability_mat = Matrix3_16d::Zero();
  child_somatic_probability = RowVector16d::Zero();
  mother_somatic_probability = RowVector16d::Zero();
  father_somatic_probability = RowVector16d::Zero();
};

bool ReadDependentData::Equals(const ReadDependentData &other) {
  bool attr_table[20] = {
    EqualsReadDataVector(read_data_vec, other.read_data_vec),
    max_elements == other.max_elements,
    sequencing_probability_mat == other.sequencing_probability_mat,
    child_somatic_probability == other.child_somatic_probability,
    mother_somatic_probability == other.mother_somatic_probability,
    father_somatic_probability == other.father_somatic_probability,
    denominator.child_zygotic_probability == other.denominator.child_zygotic_probability,
    denominator.mother_zygotic_probability == other.denominator.mother_zygotic_probability,
    denominator.father_zygotic_probability == other.denominator.father_zygotic_probability,
    denominator.child_germline_probability == other.denominator.child_germline_probability,
    denominator.parent_probability == other.denominator.parent_probability,
    denominator.root_mat == other.denominator.root_mat,
    denominator.sum == other.denominator.sum,
    numerator.child_zygotic_probability == other.numerator.child_zygotic_probability,
    numerator.mother_zygotic_probability == other.numerator.mother_zygotic_probability,
    numerator.father_zygotic_probability == other.numerator.father_zygotic_probability,
    numerator.child_germline_probability == other.numerator.child_germline_probability,
    numerator.parent_probability == other.numerator.parent_probability,
    numerator.root_mat == other.numerator.root_mat,
    numerator.sum == other.numerator.sum
  };

  if (all_of(begin(attr_table), end(attr_table), [](bool i) { return i; })) {
    return true;
  } else {
    return false;
  }
}