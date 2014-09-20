/**
 * @file bam_driver.cc
 * @author Melissa Ip
 *
 * This is the driver for parsing a bam file. The input file must contain
 * all reads for the trio and have the appropriate tags.
 *
 * SM refers to the name of the sample that identifies it as belonging to the
 * child, mother, or father. This will vary depending on the data and must be
 * known in order to parse out the read data.
 * 
 * Useful commands to view overview, header, and a certain region:
 * samtools idxstats <name>.bam
 * samtools view -H <name>.bam <chr>:<pos1>-<pos2>
 * samtools view <name>.bam <chr>:<pos1>-<pos2>
 * samtools tview -p <chr>:<pos1>-<pos2> <name>.bam
 *
 * To slice a region out and merge it with other slices:
 * samtools view -b <name>.bam <chr>:<pos1>-<pos2> > <output>.bam
 *
 * To output all header data that belong to a certain region:
 * samtools view -o /home/mip/<header>.txt -H <name>.bam <chr>:<pos1>-<pos2>
 *
 * Make read group file for bam1 to bamn, tab separated. <header>.txt can be
 * parsed by filtering out lines that begin with @RG with parse_bam_header.cc.
 *
 * Skip this step if there is only one bam file that contains the bam
 * alignments for child, mother, and father.
 * To merge separate bam files for each family member:
 * samtools merge -rh <rg>.txt <output>.bam <bam1>.bam <bamn>.bam
 *
 * To prepare the output for the index:
 * samtools sort <output>.bam <output_sorted>.bam
 *
 * To create the index:
 * samtools index <output_sorted>.bam <output>.index
 *
 * <output_sorted>.bam is accepted as input for a specified region.
 * 
 * To compile on Herschel without using cmake and include GSL and BamTools:
 * c++ -std=c++11 -L/usr/local/lib -lgsl -lgslcblas -lm -L/home/mip/novo_muta_infinite_sites_model/bamtools/lib -I/home/mip/novo_muta_infinite_sites_model/bamtools/include -lbamtools -I/home/mip/novo_muta_infinite_sites_model/bamtools/src -lbamtools-utils -I/usr/local/include -o bam_driver utility.cc read_dependent_data.cc trio_model.cc bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc em_algorithm.cc sufficient_statistics.cc bam_driver.cc
 *
 * To run this file, provide the following command line inputs:
 * ./bam_driver <output>.txt <input>.bam <child SM> <mother SM> <father SM>
 */
#include "sufficient_statistics.h"
#include "variant_visitor.h"


int main(int argc, const char *argv[]) {
/*  if (argc < 6) {
    Die("USAGE: bam_driver <output>.txt <input>.bam "
        "<child SM> <mother SM> <father SM>");
  }

  const string output_name = argv[1];
  const string input = argv[2];
  const string child_sm = argv[3];
  const string mother_sm = argv[4];
  const string father_sm = argv[5];
  const int qual_cut = 13;  // May decide to pass thresholds via command line.
  const int mapping_cut = 13;
  const double probability_cut = 0.0;  // 0.1

  BamReader reader;
  reader.Open(input);
  if (!reader.IsOpen()) {
    Die("Input file could not be opened.");
  }

  BamAlignment al;
  TrioModel params;
  PileupEngine pileup;
  SamHeader header = reader.GetHeader();
  RefVector references = reader.GetReferenceData();
  VariantVisitor *v = new VariantVisitor(references, header, params, al,
                                         output_name, child_sm, mother_sm,
                                         father_sm, qual_cut, mapping_cut,
                                         probability_cut);
  
  pileup.AddVisitor(v);
  while (reader.GetNextAlignment(al)) {
    pileup.AddAlignment(al);
  }
    
  pileup.Flush();
  reader.Close();*/

  // EM algorithm begins with initial E-Step.
  //TrioVector sites = v->sites();
  TrioVector sites;
  ReadDataVector data = {{10, 1, 0, 0}, {10, 10, 0, 0}, {40, 0, 0, 0}};
  sites.push_back(data);
  TrioModel params;
  SufficientStatistics stats(sites.size());
  stats.Update(params, sites);

  // M-Step.
  int count = 0;
  double log_likelihood = stats.log_likelihood();
  double maximized = stats.MaxSequencingErrorRate();
  while (!Equal(params.sequencing_error_rate(), maximized) &&
      count < 50) {  // Quits if converges or takes longer than 50 iteratons.
    params.set_sequencing_error_rate(maximized);  // Sets new estimate.
    stats.Clear(); // Sets to 0.
    stats.Update(params, sites);  // Loops to E-Step.
    if (stats.log_likelihood() < log_likelihood) {  // Sum of likelihood should increase and converge.
      cout << "Log likelihood is decreasing between iterations." << endl;
    }
    log_likelihood = stats.log_likelihood();
    maximized = stats.MaxSequencingErrorRate(); // Loops to M-Step.
    count++;
  }
  
  cout << "^E:\t" << params.sequencing_error_rate()
       << " after " << count << " iterations." << endl;

  return 0;
}
