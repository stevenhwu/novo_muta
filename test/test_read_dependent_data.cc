/**
 * @file test_read_dependent_data.cc
 * @author Melissa Ip
 *
 * This file tests the functions in read_dependent_data.h.
 */
#define BOOST_TEST_MODULE TestReadDependentData

#include <boost/test/unit_test.hpp>

#include "read_dependent_data.h"

BOOST_AUTO_TEST_CASE(TestGetMismatches) {
  ReadDataVector data_vec = {{20, 10, 0, 1},
                             {40, 0,  0, 0},
                             {10, 1,  1, 1}};
  ReadDependentData read_dependent_data(data_vec);
  Matrix3_16d test_matches;
  test_matches << 11.0, 1.0,  11.0, 10.0,
                  1.0,  21.0, 21.0, 20.0,
                  11.0, 21.0, 31.0, 30.0,
                  10.0, 20.0, 30.0, 30.0,

                  0.0,  0.0,  0.0,  0.0,
                  0.0,  40.0, 40.0, 40.0,
                  0.0,  40.0, 40.0, 40.0,
                  0.0,  40.0, 40.0, 40.0,

                  3.0,  2.0,  2.0,  2.0,
                  2.0,  12.0, 11.0, 11.0,
                  2.0,  11.0, 12.0, 11.0,
                  2.0,  11.0, 11.0, 12.0;

  Matrix3_16d matches = read_dependent_data.GetMismatches();
  BOOST_CHECK(matches == test_matches);
}

BOOST_AUTO_TEST_CASE(TestGetHeterozygousMatches) {
  ReadDataVector data_vec = {{20, 10, 0, 1},
                             {40, 0,  0, 0},
                             {10, 1,  1, 1}};
  ReadDependentData read_dependent_data(data_vec);
  Matrix3_16d test_matches;
  test_matches << 0.0,  30.0, 20.0, 21.0,
                   30.0, 0.0,  10.0, 11.0,
                   20.0, 10.0, 0.0,  1.0,
                   21.0, 11.0, 1.0,  0.0,

                   0.0,  40.0, 40.0, 40.0,
                   40.0, 0.0,  0.0,  0.0,
                   40.0, 0.0,  0.0,  0.0,
                   40.0,  0.0,  0.0,  0.0,

                   0.0,  11.0, 11.0, 11.0,
                   11.0, 0.0,  2.0,  2.0,
                   11.0, 2.0,  0.0,  2.0,
                   11.0, 2.0,  2.0,  0.0;

  Matrix3_16d matches = read_dependent_data.GetHeterozygousMatches();
  BOOST_CHECK(matches == test_matches);
}

BOOST_AUTO_TEST_CASE(TestGetHomozygousMatches) {
  ReadDataVector data_vec = {{20, 10, 0, 1},
                             {40, 0,  0, 0},
                             {10, 1,  1, 1}};
  ReadDependentData read_dependent_data(data_vec);
  Matrix3_16d test_matches;
  test_matches << 20.0, 0.0,  0.0, 0.0,
                  0.0,  10.0, 0.0, 0.0,
                  0.0,  0.0,  0.0, 0.0,
                  0.0,  0.0,  0.0, 1.0,

                  40.0, 0.0,  0.0, 0.0,
                  0.0,  0.0,  0.0, 0.0,
                  0.0,  0.0,  0.0, 0.0,
                  0.0,  0.0,  0.0, 0.0,

                  10.0, 0.0,  0.0, 0.0,
                  0.0,  1.0,  0.0, 0.0,
                  0.0,  0.0,  1.0, 0.0,
                  0.0,  0.0,  0.0, 1.0;

  Matrix3_16d matches = read_dependent_data.GetHomozygousMatches();
  BOOST_CHECK(matches == test_matches);
}