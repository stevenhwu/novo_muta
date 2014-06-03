/**
 * @file bam_driver.cc
 * @author Melissa Ip
 *
 * This is a test driver for parsing a bam file.
 *
 * To compile on Herschel and include GSL and BamTools:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -L/home/mip/novo_muta_infinite_sites_model/bamtools/lib -I/home/mip/novo_muta_infinite_sites_model/bamtools/include -lbamtools -lbamtools-utils -I/usr/local/include -o bam_driver utility.cc bam_driver.cc
 */
#include "api/BamReader.h"

#include "utility.h"

using namespace BamTools;


int main(int argc, const char *argv[]) {
  if (argc < 2) {
    Die("USAGE: bam_driver <input>.bam");
  }

  const string file_name = argv[1];

  BamReader fin;
  fin.Open(file_name);

  if (!fin.IsOpen()) {
    Die("Input file could not be opened.");
  }

  BamAlignment al;
  while (fin.GetNextAlignment(al)) {
    cout << al.QueryBases << endl;
  }

  fin.Close();

  return 0;
}
