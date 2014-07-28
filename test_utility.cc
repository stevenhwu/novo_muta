/**
 * @file test_utility.cc
 * @author Melissa Ip
 *
 * This file tests the functions in utility.h.
 *
 * To compile on Herschel, go to the build directory and call:
 * cmake ..
 * make
 */
#define BOOST_TEST_MODULE TestUtility

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "utility.h"


// BOOST_AUTO_TEST_CASE(TestGenotypeNumIndex) {
//   Matrix16_2i mat = GenotypeNumIndex();
//   for (int i = 0; i < kNucleotideCount; ++i) {
//     BOOST_CHECK(mat(i, 0) == 0);
//     BOOST_CHECK(mat(i+4, 0) == 1);
//     BOOST_CHECK(mat(i+8, 0) == 2);
//     BOOST_CHECK(mat(i+12, 0) == 3);
//     BOOST_CHECK(mat(i, 1) == i);
//     BOOST_CHECK(mat(i+4, 1) == i);
//     BOOST_CHECK(mat(i+8, 1) == i);
//     BOOST_CHECK(mat(i+12, 1) == i);
//   }
// }

BOOST_AUTO_TEST_CASE(TestZeroMatrix16_16_4d) {
  Matrix16_16_4d mat = ZeroMatrix16_16_4d();
  BOOST_CHECK(mat.rows() == kGenotypeCount);
  BOOST_CHECK(mat.cols() == kGenotypeCount);
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      BOOST_CHECK(mat(i, j).size() == kNucleotideCount);
      for (int k = 0; k < kNucleotideCount; ++k) {
        BOOST_CHECK(mat(i, j)(k) == 0.0);
      }
    }
  }
}

// BOOST_AUTO_TEST_CASE(TestPrintMatrix16_16_4d) {}

BOOST_AUTO_TEST_CASE(TestTwoParentCounts) {
  Matrix16_16_4d mat = TwoParentCounts();
  RowVector4d vec = RowVector4d::Zero();
  for (int i = 0; i < kGenotypeCount; ++i) {
    for (int j = 0; j < kGenotypeCount; ++j) {
      vec(i / kNucleotideCount)++;
      vec(i % kNucleotideCount)++;
      vec(j / kNucleotideCount)++;
      vec(j % kNucleotideCount)++;
      for (int k = 0; k < kNucleotideCount; ++k) {
        RowVector4d nucleotide_counts = mat(i, j);
        BOOST_CHECK(nucleotide_counts(k) == vec(k));
      }
      vec = RowVector4d::Zero();
    }
  }
}

BOOST_AUTO_TEST_CASE(TestEqualsReadData) {
  ReadData a = {1, 1, 1, 1};
  ReadData b = {1, 1, 1, 1};
  ReadData c = {2, 2, 2, 2};
  BOOST_CHECK(EqualsReadData(a, b));
  BOOST_CHECK(!EqualsReadData(a, c));
}

BOOST_AUTO_TEST_CASE(TestEqualsReadDataVector) {
  ReadDataVector a = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
  ReadDataVector b = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
  ReadDataVector c = {{2, 2, 2, 2}, {2, 2, 2, 2}, {2, 2, 2, 2}};
  BOOST_CHECK(EqualsReadDataVector(a, b));
  BOOST_CHECK(!EqualsReadDataVector(a, c));
}

// BOOST_AUTO_TEST_CASE(TestPrintReadData) {}

// BOOST_AUTO_TEST_CASE(TestPrintReadDataVector) {}

BOOST_AUTO_TEST_CASE(TestEnumerateNucleotideCounts) {
  ReadDataVector vec = EnumerateNucleotideCounts(kNucleotideCount);
  BOOST_CHECK(vec.size() == 624);
}

BOOST_AUTO_TEST_CASE(TestGetTrioVector) {
  TrioVector trio_vec = GetTrioVector(kNucleotideCount);
  BOOST_CHECK(trio_vec.size() == kTrioCount);
  ReadDataVector vec = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  BOOST_CHECK(IndexOfReadDataVector(vec, trio_vec) == -1);
}

BOOST_AUTO_TEST_CASE(TestIndexOfReadDataVector) {
  ReadDataVector a = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  ReadDataVector b = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
  ReadDataVector c = {{2, 2, 2, 2}, {2, 2, 2, 2}, {2, 2, 2, 2}};
  TrioVector trio_vec = {a, b};
  BOOST_CHECK(IndexOfReadDataVector(b, trio_vec) == 1);
  BOOST_CHECK(IndexOfReadDataVector(c, trio_vec) == -1);
}

BOOST_AUTO_TEST_CASE(TestIsInVector) {
  RowVector4d vec;
  vec << 0.0, 1.0, 2.0, 3.0;
  BOOST_CHECK(IsInVector(vec, 1.0));
  BOOST_CHECK(!IsInVector(vec, 4.0));
}

BOOST_AUTO_TEST_CASE(TestIsAlleleInParentGenotype) {
  for (int i = 0; i < kGenotypeCount; ++i) {
    BOOST_CHECK(IsAlleleInParentGenotype(i / kNucleotideCount, i));
    BOOST_CHECK(IsAlleleInParentGenotype(i % kNucleotideCount, i));
  }
}

BOOST_AUTO_TEST_CASE(TestDirichletMultinomialLog) {}

BOOST_AUTO_TEST_CASE(TestKroneckerProduct) {
  Matrix4_16d mat;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 16; ++j) {
      mat(i, j) = j;
    }
  }
  Matrix16_256d kron = KroneckerProduct(mat);
  BOOST_REQUIRE(kron.rows() == 16);
  BOOST_REQUIRE(kron.cols() == 256);
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 256; ++j) {
      BOOST_CHECK(kron(i, j) == mat(i/4, j/16) * mat(i%4, j%16));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestKroneckerProduct2) {
  Matrix4d mat;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mat(i, j) = j;
    }
  }
  Matrix16_16d kron = KroneckerProduct(mat);
  BOOST_REQUIRE(kron.rows() == 16);
  BOOST_REQUIRE(kron.cols() == 16);
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {
      BOOST_CHECK(kron(i, j) == mat(i/4, j/4) * mat(i%4, j%4));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestKroneckerProduct3) {
  RowVector16d vec;
  for (int i = 0; i < 16; ++i) {
    vec(i) = i;
  }
  RowVector256d kron = KroneckerProduct(vec, vec);
  BOOST_REQUIRE(kron.rows() == 1);
  BOOST_REQUIRE(kron.cols() == 256);
  for (int i = 0; i < 256; ++i) {
    BOOST_CHECK(kron(i) == vec(i/16) * vec(i%16));
  }
}

BOOST_AUTO_TEST_CASE(TestEqual) {
  BOOST_CHECK(Equal(0.0, 0.0));
  BOOST_CHECK(!Equal(0.0, 1.0));
  BOOST_CHECK(Equal(0.0, 2.22e-16));
  BOOST_CHECK(!Equal(0.0, 2.23e-16));
}

// BOOST_AUTO_TEST_CASE(TestDie) {}