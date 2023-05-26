/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <filesystem>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
using std::string;
#include <atomic>
#include <map>
// namespace fs = std::filesystem;

#include <spdlog/spdlog.h>

// Alternative interval tree implemented base on red-black tree.
#include "ygg.hpp"

#include "bamUtils.h"
#include "geneBuilder.h"

using MyInterval = std::pair< int, int >;
template < class Node > class NodeTraits : public ygg::ITreeNodeTraits< Node >
{
public:
    using key_type = int;
    static int get_lower(const Node& node)
    {
        return node.lower;
    }
    static int get_upper(const Node& node)
    {
        return node.upper;
    }
    static int get_lower(const MyInterval& i)
    {
        return std::get< 0 >(i);
    }
    static int get_upper(const MyInterval& i)
    {
        return std::get< 1 >(i);
    }
};

class Node : public ygg::ITreeNodeBase< Node, NodeTraits< Node > >
{
public:
    int         upper;
    int         lower;
    GeneFromGTF value;
};

using MyTree = ygg::IntervalTree< Node, NodeTraits< Node > >;

enum LocusFunction
{
    NONE,
    INTERGENIC,
    INTRONIC,
    CODING
};

class TagReadsWithGeneExon
{
public:
    TagReadsWithGeneExon(string annotation_filename)
        : annotation_filename(annotation_filename), multiMap(false)
    {
        total_reads         = 0;
        map_reads           = 0;
        antisense_reads     = 0;
        exonic_reads        = 0;
        intronic_reads      = 0;
        intergenic_reads    = 0;
        transcriptome_reads = 0;
    }
    ~TagReadsWithGeneExon() {}

    int makeOverlapDetector();

    // record, contig: input parameters
    // gene: output genename
    // strand: output strand, 1 => positive 2 => negative
    // locus: output locus function, 1 => intergenic 2 => intronic 3 => exonic
    int setAnnotation(BamRecord& record, const std::string& contig, string& gene,
                      unsigned int* strand, unsigned int* locus);

    // get exon length and intron length, for multi-mapping
    int setAnnotation(BamRecord& record, const std::string& contig, int& overlap,
                      int& locus);

    std::string dumpMetrics();

    std::unordered_map< std::string, unsigned > getGeneID()
    {
        return geneID;
    }
    std::vector< std::string > getGeneName()
    {
        return geneName;
    }
    std::map< std::string, std::vector< pair< std::string, unsigned int > > > getGeneEnd()
    {
        return geneEnd;
    }

    void setMultiMap()
    {
        multiMap = true;
    }

private:
    // Calculate locus function by comparing then overlap length of alignment
    // block and exons/introns
    std::pair< int, LocusFunction >
    getLocusFunction(std::vector< const GeneFromGTF* >& result,
                     std::vector< AlignmentBlock >& alignmentBlock, int& overlap);
    std::pair< int, LocusFunction >
    getLocusFunction(std::vector< const GeneFromGTF* >& result,
                     std::vector< AlignmentBlock >&     alignmentBlock);
    // Set annotation by overlap length
    int annoByOverlapLen(BamRecord& record, const std::string& contig, string& gene,
                         unsigned int* strand, unsigned int* locus);

private:
    string annotation_filename;

    std::unordered_map< std::string, MyTree > mytrees;
    std::vector< Node >                       nodes;

    // Annotation statics
    std::atomic< size_t > total_reads;
    std::atomic< size_t > map_reads;
    std::atomic< size_t > antisense_reads;
    std::atomic< size_t > exonic_reads;
    std::atomic< size_t > intronic_reads;
    std::atomic< size_t > intergenic_reads;
    std::atomic< size_t > transcriptome_reads;

    std::unordered_map< std::string, unsigned >                               geneID;
    std::vector< std::string >                                                geneName;
    std::map< std::string, std::vector< pair< std::string, unsigned int > > > geneEnd;

    bool multiMap;
};
