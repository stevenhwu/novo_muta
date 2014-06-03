/**
 * @file bam_driver.cc
 * @author Melissa Ip
 *
 * This is a test driver for parsing a bam file.
 *
 * To compile on Herschel and include GSL and BamTools:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -L/home/mip/novo_muta_infinite_sites_model/bamtools/lib -I/home/mip/novo_muta_infinite_sites_model/bamtools/include -lbamtools -I/home/mip/novo_muta_infinite_sites_model/bamtools/src -lbamtools-utils -I/usr/local/include -o bam_driver utility.cc read_dependent_data.cc trio_model.cc bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc bam_driver.cc
 */
#include "variant_visitor.h"


int main(int argc, const char *argv[]) {
  if (argc < 2) {
    Die("USAGE: bam_driver <input>.bam");
  }

  const string file_name = argv[1];

  BamReader reader;
  reader.Open(file_name);

  if (!reader.IsOpen()) {
    Die("Input file could not be opened.");
  }

  RefVector references = reader.GetReferenceData();
  SamHeader header = reader.GetHeader();
  PileupEngine pileup;
  TrioModel params;
  BamAlignment al;
  VariantVisitor *v = new VariantVisitor(references, header, params, al, 13, 13, 0.1);
  pileup.AddVisitor(v);
   
  while (reader.GetNextAlignment(al)) {
    pileup.AddAlignment(al);
  }
    
  pileup.Flush();
  reader.Close();

  return 0;
}
