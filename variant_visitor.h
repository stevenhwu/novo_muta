/**
 * @file variant_visitor.h
 * @author Melissa Ip
 *
 * The VariantVisitor class inherits the PileupVisitor class in BamTools.
 * It is used to parse bam files and retrieve sequencing reads, which are
 * accepted by the trio model.
 *
 * This class is referenced from:
 * https://github.com/dwinter/accuMUlate
 */
#ifndef VARIANT_VISITOR_H
#define VARIANT_VISITOR_H

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"

#include "trio_model.h"

using namespace BamTools;


/**
 * VariantVisitor class header. See top of file for a complete description.
 */
class VariantVisitor : public PileupVisitor {
 public:
  VariantVisitor(const RefVector &references,
                 const SamHeader &header,
                 const TrioModel &params,
                 BamAlignment &al,
                 string output_name,  // Unused.
                 string child_sm,
                 string mother_sm,
                 string father_sm,
                 int qual_cut,
                 int mapping_cut,
                 double probability_cut);
  ~VariantVisitor() { }
  void Visit(const PileupPosition &pileupData);
  TrioVector sites() const;

 private:
  RefVector references_;
  SamHeader header_;
  TrioModel params_;
  BamAlignment &al_;
  string output_name_;
  string child_sm_;
  string mother_sm_;
  string father_sm_;
  int qual_cut_;
  int mapping_cut_;
  double probability_cut_;
  TrioVector sites_;
};

#endif
