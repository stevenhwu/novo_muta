/**
 * @file test_pileup_utility.cc
 * @author Melissa Ip
 *
 * This file tests the functions in pileup_utility.h.
 *
 * To compile on Herschel, go to the build directory and call:
 * cmake ..
 * make
 */
#define BOOST_TEST_MODULE TestPileupUtility

#include <boost/test/unit_test.hpp>

#include "pileup_utility.h"


// BOOST_AUTO_TEST_CASE(TestTrimHeader) {}

BOOST_AUTO_TEST_CASE(TestGetSequence) {
  string line = "22\t34907918\tG\t83\t"
                ".TT,T.TT.tt,..t..TT..tT.tt,..tt.t.TttT,Tt,tT,tt,T,T,Tt,t.TtTTTTTt,.,.TTT.,,T..ttT.,\t"
                "DEBB/>EFBFFC<EFEEFFEDFBCFECEEGEEFEFBGE>FGCGFC;GCFC;CEGCFEEGEEEFEGCECEEFED;C>CAFFDC@\n";
  BOOST_CHECK(GetSequence(line).find('\n') == string::npos);
  string invalid_line = "21\t9411191\tN\t1\t^>T C\n";
  BOOST_CHECK(GetSequence(invalid_line).empty());
}

BOOST_AUTO_TEST_CASE(TestGetReadData) {
  string line = "22\t34907918\tG\t83\t"
                ".TT,T.TT.tt,..t..TT..tT.tt,..tt.t.TttT,Tt,tT,tt,T,T,Tt,t.TtTTTTTt,.,.TTT.,,T..ttT.,\t"
                "DEBB/>EFBFFC<EFEEFFEDFBCFECEEGEEFEFBGE>FGCGFC;GCFC;CEGCFEEGEEEFEGCECEEFED;C>CAFFDC@";
  ReadData data = GetReadData(line);
  ReadData test_data = {0, 0, 36, 47};
  BOOST_CHECK(EqualsReadData(test_data, data));
}

BOOST_AUTO_TEST_CASE(TestGetProbability) {
  TrioModel params;
  string child_line = "22\t34907918\tG\t83\t"
                      ".TT,T.TT.tt,..t..TT..tT.tt,..tt.t.TttT,Tt,tT,tt,T,T,Tt,t.TtTTTTTt,.,.TTT.,,T..ttT.,\t"
                      "DEBB/>EFBFFC<EFEEFFEDFBCFECEEGEEFEFBGE>FGCGFC;GCFC;CEGCFEEGEEEFEGCECEEFED;C>CAFFDC@";
  string mother_line = "22\t34907918\tG\t60\t"
                       ".$,,,,.,.,....,..,,.,,.,.,..,,,,..,.,.....,.,.....,......,,.,\t"
                       "@@BABEBCBBCDDCEECBDBC@CDCDDBBDCABBDCDCD@DBDBCCD1DCCC@CCCCBBB";
  string father_line = "22\t34907918\tG\t68\t"
                       ".....,.,,.,,....,,,,,.,,,,,.,..,,,.,,..,,.,,,..,,,.,.,..,.,,,..,.,.,\t"
                       "DADADBACADBBDEEDCCBCBECCCCBEBDDCCBCCC9DBBEBCCEDCCDEDDDC<BCBBBDCBB@?5";
  double probability = GetProbability(params, child_line, mother_line, father_line);
  BOOST_CHECK(probability >= 0.0 && probability <= 1.0);
}

// BOOST_AUTO_TEST_CASE(TestProcessPileup) {}