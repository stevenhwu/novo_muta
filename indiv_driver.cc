/**
 * @file indiv_driver.cc
 * @author Melissa Ip
 *
 * This allows the user to pass in an individual site that is converted to a
 * ReadDataVector. Using the trio model, the probability of mutation is
 * calculated and printed.
 */
#include <sstream>

#include "trio_model.h"


int main() {
  string child_str;
  string mother_str;
  string father_str;
  uint16_t A = 0;
  uint16_t C = 0;
  uint16_t G = 0;
  uint16_t T = 0;

  while (true) {
    cout << "Please enter the child sequencing read <A C G T>. Example: 0 1 2 3" << endl;
    getline(cin, child_str);
    stringstream child(child_str);
    child >> A;
    child >> C;
    child >> G;
    child >> T;
    ReadData child_read = {A, C, G, T};

    cout << endl << "Please enter the mother sequencing read <A,C,G,T>. Example: 0,1,2,3" << endl;
    getline(cin, mother_str);
    stringstream mother(mother_str);
    mother >> A;
    mother >> C;
    mother >> G;
    mother >> T;
    ReadData mother_read = {A, C, G, T};

    cout << endl << "Please enter the father sequencing read <A,C,G,T>. Example: 0,1,2,3" << endl;
    getline(cin, father_str);
    stringstream father(father_str);
    father >> A;
    father >> C;
    father >> G;
    father >> T;
    ReadData father_read = {A, C, G, T};

    TrioModel params;
    ReadDataVector data = {child_read, mother_read, father_read};
    cout << endl << "P(Mut):\t" << params.MutationProbability(data) << endl << endl;
  }

  return 0;
}
