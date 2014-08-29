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
 * nucleotide bases from read groups (RG). Creates a ReadDataVector.
 *
 * @param  pileupData  Pileup data.
 */
void VariantVisitor::Visit(const PileupPosition &pileupData) {
  ReadDataVector data_vec = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  string chr = references_[pileupData.RefId].RefName;
  uint64_t pos = pileupData.Position;
  string tag_id;
  bool is_qual = false;

  for (auto it = begin(pileupData.PileupAlignments);
      it != end(pileupData.PileupAlignments); ++it) {
    if (it->Alignment.MapQuality >= mapping_cut_) {
      const int *pos = &it->PositionInAlignment;
      uint16_t bqual = static_cast<short>(it->Alignment.Qualities[*pos]) - 33;
      if (bqual >= qual_cut_) {
        string sm;
        int i = 0;
        for (SamReadGroupConstIterator rg_it = header_.ReadGroups.ConstBegin();
            rg_it != header_.ReadGroups.ConstEnd(); ++rg_it) {
          sm = rg_it->Sample;
          if (sm.compare(child_sm_) == 0) {
            i = 0;  // Hard coded as implemented in TrioModel.
          } else if (sm.compare(mother_sm_) == 0) {
            i = 1;
          } else if (sm.compare(father_sm_) == 0) {
            i = 2;
          } else {
            cout << "ERROR: " << sm << " does not match." << endl;
          }
        
          char base = it->Alignment.QueryBases[*pos];
          uint16_t base_idx = ToNucleotideIndex(base);
          if (base_idx >= 0 && base_idx < 4) {
            is_qual = true;
            data_vec[i].reads[base_idx]++;
          }
        }
      }
    }
  }

  if (is_qual) {
    sites_.push_back(data_vec);
  }
}

TrioVector VariantVisitor::sites() const {
  return sites_;
}