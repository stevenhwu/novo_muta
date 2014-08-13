/**
 * @file test_sufficient_statistics.cc
 * @author Melissa Ip
 *
 * This file tests the functions in sufficient_statistics.h.
 *
 * To compile on Herschel, go to the build directory and call:
 * cmake ..
 * make
 */
#define BOOST_TEST_MODULE TestSufficientStatistics

#include <boost/test/unit_test.hpp>

#include "sufficient_statistics.h"


BOOST_AUTO_TEST_CASE(TestSufficientStatistics) {
  SufficientStatistics stats(1.0);
  BOOST_CHECK(stats.e() == 0.0);
  BOOST_CHECK(stats.hom() == 0.0);
  BOOST_CHECK(stats.het() == 0.0);
  BOOST_CHECK(stats.som() == 0.0);
  BOOST_CHECK(stats.germ() == 0.0);
  BOOST_CHECK(stats.n_s() == 1.0);
}

BOOST_AUTO_TEST_CASE(TestMaxGermlineMutationRate) {
  SufficientStatistics stats(1.0);
  double s_germ = stats.MaxGermlineMutationRate();
  BOOST_CHECK(s_germ == 0.0);
}

BOOST_AUTO_TEST_CASE(TestMaxSomaticMutationRate) {
  SufficientStatistics stats(1.0);
  double s_som = stats.MaxSomaticMutationRate();
  BOOST_CHECK(s_som == 0.0);
}

BOOST_AUTO_TEST_CASE(TestMaxSequencingErrorRate) {
  SufficientStatistics stats(1.0);
  stats.set_hom(1.0);
  stats.set_het(1.0);
  double s_e = stats.MaxSequencingErrorRate();
  BOOST_CHECK(s_e == 0.0);
}

BOOST_AUTO_TEST_CASE(TestUpdate) {
  TrioModel params;
  SufficientStatistics stats(1.0);
  TrioVector vec;
  ReadDataVector data = {{40, 0, 0, 0},
                         {40, 0, 0, 0},
                         {40, 0, 0, 0}};
  vec.push_back(data);
  stats.Update(params, vec);
  BOOST_CHECK(stats.som() > 0.0);
  BOOST_CHECK(stats.germ() > 0.0);
  BOOST_CHECK(stats.e() > 0.0);
  BOOST_CHECK(stats.hom() > 0.0);
  BOOST_CHECK(stats.het() > 0.0);
}

BOOST_AUTO_TEST_CASE(TestClear) {
  SufficientStatistics stats(1.0);
  stats.set_e(1.0);
  stats.set_hom(1.0);
  stats.set_het(1.0);
  stats.set_som(1.0);
  stats.set_germ(1.0);
  stats.Clear();
  BOOST_CHECK(stats.e() == 0.0);
  BOOST_CHECK(stats.hom() == 0.0);
  BOOST_CHECK(stats.het() == 0.0);
  BOOST_CHECK(stats.som() == 0.0);
  BOOST_CHECK(stats.germ() == 0.0);
}

BOOST_AUTO_TEST_CASE(TestIsNan) {
  SufficientStatistics stats(1.0);
  BOOST_CHECK(!stats.IsNan());
  stats.set_het(0.0 / 0.0);
  BOOST_CHECK(stats.IsNan());
}

// BOOST_AUTO_TEST_CASE(TestPrint) {}