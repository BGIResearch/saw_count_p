/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "geneFromGTF.h"
#include "gtfRecord.h"

typedef std::unordered_map< std::string, std::vector< GTFRecord > >::const_iterator
                                                                    GTFIterator;
typedef std::unordered_map< std::string, std::vector< GTFRecord > > GTFMapForTranscript;

class GeneBuilder
{
public:
    GeneBuilder() {}
    // Main entrance, given iterator of gtf records, return one or more gene
    // structure
    std::vector< GeneFromGTF > makeGene(GTFIterator& iter);

private:
    // Check if the genes in different contigs
    bool geneInDiffChrs(std::vector< GTFRecord >& gtfRecords);
    // Grouping the records with transcript id
    GTFMapForTranscript gatherByTransciptID(std::vector< GTFRecord >& gtfRecords);
    // There is not exists different contigs, so do the work
    GeneFromGTF makeGeneWithTranscripts(std::vector< GTFRecord >& gtfRecords);
    // Add gene
    GeneFromGTF makeGene(std::vector< GTFRecord >& gtfRecords);
    // Add transcript
    void addTranscriptToGene(GeneFromGTF& gene, std::vector< GTFRecord >& gtfRecords);

    void validateGTFRecord(GTFRecord& record, GeneFromGTF& gene);
};
