/**
 * @file variant_visitor.cc
 * @author Melissa Ip
 *
 * This file contains the implementation of the VariantVisitor class.
 *
 * See top of variant_visitor.h for a complete description.
 */
#include "variant_visitor.h"


/**
 * Constructor calls base constructor PileupVisitor.
 */
VariantVisitor::VariantVisitor(const RefVector& references,
              								 const SamHeader& header,
								               const TrioModel& params,
								               BamAlignment& al,
								               int qual_cut,
								               int mapping_cut,
								               double prob_cut) :
    PileupVisitor(),
    references_{references}, header_{header}, al_{al}, qual_cut_{qual_cut},
    prob_cut_{prob_cut}, mapping_cut_{mapping_cut} {
}

/**
 * Visits a bam alignment at a certain position and counts the number of
 * nucleotide bases. Creates a ReadDataVector and calculates the probability of
 * mutation.
 *
 * @param  pileupData  Pileup data.
 */
void VariantVisitor::Visit(const PileupPosition& pileupData) {
	ReadDataVector data_vec = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
	string chr = references_[pileupData.RefId].RefName;
  uint64_t pos = pileupData.Position;
  string tag_id;
  for (auto it = begin(pileupData.PileupAlignments);
  		it != end(pileupData.PileupAlignments); ++it) {
    if (it->Alignment.MapQuality >= mapping_cut_) {
      const int *pos = &it->PositionInAlignment;
      uint16_t bqual = static_cast<short>(it->Alignment.Qualities[*pos]) - 33;
      if (bqual >= qual_cut_) {
        it->Alignment.GetTag("RG", tag_id);
        string sm = header_.ReadGroups[tag_id].Sample;
        cout << sm << endl;

        char base = it->Alignment.QueryBases[*pos];
        uint16_t bindex = VariantVisitor::ToNucleotideIndex(base);
        if (bindex < 4) {
          //data_vec[/*compare sm to family member].reads[bindex]++;
        }
	    }
    }
	}

  // double probability = params_.MutationProbability(data_vec);
  // if (probability >= prob_cut_) {
  //   cout << chr << '\t' << pos << '\t' << probability << endl;
  // }
}

/**
 * Returns numeric nucleotide index given a char.
 *
 * @param  b Char.
 * @return   Numeric nucleotide index.
 */
uint16_t VariantVisitor::ToNucleotideIndex(char b) {
  switch(b) {        
  case 'A':
  case 'a':    
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  case '-':
  case 'N':
    return -1 ;
  default:  // Unknown base.
	  return -1;
  }
}