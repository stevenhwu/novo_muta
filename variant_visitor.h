/**
 * @file variant_visitor.h
 * @author Melissa Ip
 *
 * The VariantVisitor class inherits the PileupVisitor class in BamTools.
 * It is used to parse bam files and retrieve sequencing reads, which are
 * accepted by the trio model.
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
  VariantVisitor(const RefVector& references,
                 const SamHeader& header,
                 const TrioModel& params,
                 BamAlignment& al,
                 int qual_cut,
                 int mapping_cut,
                 double prob_cut);
  ~VariantVisitor() { }
  void Visit(const PileupPosition& pileupData);
  uint16_t ToNucleotideIndex(char b);

 private:
  RefVector references_;
  SamHeader header_;
  BamAlignment& al_;
  TrioModel params_;
  int qual_cut_;
  int mapping_cut_;
  double prob_cut_;
};

#endif
