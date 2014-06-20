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
VariantVisitor::VariantVisitor(const RefVector &references,
              								 const SamHeader &header,
								               const TrioModel &params,
								               BamAlignment &al,
                               string output_name,
                               string child_sm,
                               string mother_sm,
                               string father_sm,
								               int qual_cut,
								               int mapping_cut,
								               double probability_cut)
    : PileupVisitor(),
      references_{references}, header_{header}, al_{al},
      output_name_{output_name}, child_sm_{child_sm}, mother_sm_{mother_sm},
      father_sm_{father_sm}, qual_cut_{qual_cut}, mapping_cut_{mapping_cut},
      probability_cut_{probability_cut} {
}

/**
 * Visits a bam alignment at a certain position and counts the number of
 * nucleotide bases from read groups (RG). Creates a ReadDataVector and
 * calculates the probability of mutation, which is printed to a text file.
 *
 * @param  pileupData  Pileup data.
 */
void VariantVisitor::Visit(const PileupPosition &pileupData) {
  ofstream fout(output_name_);
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

        int i = 0;
        if (sm.compare(child_sm_) == 0) {
          i = 0;  // Hard coded as implemented in TrioModel.
        } else if (sm.compare(mother_sm_) == 0) {
          i = 1;
        } else if (sm.compare(father_sm_) == 0) {
          i = 2;
        }

        char base = it->Alignment.QueryBases[*pos];
        uint16_t base_idx = VariantVisitor::ToNucleotideIndex(base);
        if (base_idx >= 0 && base_idx < 4) {
          data_vec[i].reads[base_idx]++;
        }
	    }
    }
	}
  
  double probability = params_.MutationProbability(data_vec);
  if (probability >= probability_cut_) {
    fout << chr << '\t' << pos << '\t' << probability << endl;
  }
  fout.close();
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