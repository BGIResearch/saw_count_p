/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include <algorithm>
#include <fstream>
#include <sstream>

#include "tagReadsWithGeneExon.h"

#include "annotationException.h"
#include "bamUtils.h"
#include "gtfReader.h"
#include <libx/Timer.hpp>

struct CmpLocusExon
{
    bool operator()(const pair< int, int >& a, const pair< int, int >& b)
    {
        return a.first < b.first || (a.first == b.first && a.second < b.second);
    }
};

struct CmpLocusIntron
{
    bool operator()(const pair< int, int >& a, const pair< int, int >& b)
    {
        return a.second < b.second;
    }
};

std::pair< int, LocusFunction >
TagReadsWithGeneExon::getLocusFunction(std::vector< const GeneFromGTF* >& result,
                                       std::vector< AlignmentBlock >&     alignmentBlock)
{
    // Calculate LocusFunction for each overlapped Gene
    int gidx = 0;

    vector< LocusFunction >    locusResult(result.size());
    vector< pair< int, int > > total_cnts(result.size());
    for (unsigned j = 0; j < result.size(); ++j)
    {
        const GeneFromGTF* gene     = result[j];
        pair< int, int >   max_cnts = { 0, 0 };  // <exon numbers, intron numbers>
        // For every block, find the confidently result
        for (auto& b : alignmentBlock)
        {
            pair< int, int > cnts = { 0, 0 };
            for (auto& p : gene->getTranscripts())
            {
                p.second.assignLocusFunction(b.getReferenceStart(), b.getLength(), cnts);
            }
            max_cnts.first += cnts.first;
            max_cnts.second += cnts.second;
        }  // end Block

        int len = 0;
        for (auto& b : alignmentBlock)
            len += b.getLength();
        LocusFunction locusFunction(LocusFunction::INTERGENIC);
        if (max_cnts.first * 2 >= len)
            locusFunction = LocusFunction::CODING;
        else if (max_cnts.second * 2 >= len)
            locusFunction = LocusFunction::INTRONIC;

        locusResult[j] = locusFunction;
        total_cnts[j]  = max_cnts;
        // cout<<j<<" "<<max_cnts.first<<" "<<max_cnts.second<<
        //     " "<<len<<" "<<locusFunction<<endl;
    }  // end Gene

    // Pick the most confidently
    if (result.size() > 1)
    {
        std::vector< int > genes;
        for (unsigned i = 0; i < result.size(); ++i)
        {
            if (genes.empty())
                genes.push_back(i);
            else
            {
                if (locusResult[i] > locusResult[genes[0]])
                {
                    // Find better gene
                    genes.clear();
                    genes.push_back(i);
                }
                else if (locusResult[i] == locusResult[genes[0]])
                {
                    genes.push_back(i);
                }
            }
        }

        std::unordered_map< int, LocusFunction > tmp_map;
        if (genes.size() == 1)
        {
            // Lucky, only one gene
            gidx = genes.front();
        }
        else
        {
            // Pick one gene using the overlap numbers of exon and intro
            vector< pair< int, int > > filter_cnts;
            for (auto& pos : genes)
            {
                filter_cnts.push_back(total_cnts[pos]);
            }
            if (locusResult[genes.front()] == LocusFunction::CODING)
            {
                auto iter = std::max_element(filter_cnts.begin(), filter_cnts.end(),
                                             CmpLocusExon());
                gidx      = genes[iter - filter_cnts.begin()];
            }
            else if (locusResult[genes.front()] == LocusFunction::INTRONIC)
            {
                auto iter = std::max_element(filter_cnts.begin(), filter_cnts.end(),
                                             CmpLocusIntron());
                gidx      = genes[iter - filter_cnts.begin()];
            }
        }
    }

    return { gidx, locusResult[gidx] };
}

std::pair< int, LocusFunction >
TagReadsWithGeneExon::getLocusFunction(std::vector< const GeneFromGTF* >& result,
                                       std::vector< AlignmentBlock >&     alignmentBlock,
                                       int&                               overlap)
{
    // Calculate LocusFunction for each overlapped Gene
    int gidx = 0;

    vector< LocusFunction > locusResult(result.size());
    vector< int >           total_cnts(result.size());
    for (unsigned j = 0; j < result.size(); ++j)
    {
        const GeneFromGTF* gene     = result[j];
        pair< int, int >   max_cnts = { 0, 0 };  // <exon numbers, intron numbers>
        // For every block, find the confidently result
        for (auto& b : alignmentBlock)
        {
            pair< int, int > cnts = { 0, 0 };
            for (auto& p : gene->getTranscripts())
            {
                p.second.assignLocusFunction(b.getReferenceStart(), b.getLength(), cnts);
            }
            max_cnts.first += cnts.first;
            max_cnts.second += cnts.second;
        }  // end Block

        int len = 0;
        for (auto& b : alignmentBlock)
            len += b.getLength();
        LocusFunction locusFunction(LocusFunction::INTERGENIC);
        int           overlapLen = 0;
        if (max_cnts.first * 2 >= len)
        {
            locusFunction = LocusFunction::CODING;
            overlapLen    = max_cnts.first;
        }
        else if (max_cnts.second * 2 >= len)
        {
            locusFunction = LocusFunction::INTRONIC;
            overlapLen    = max_cnts.second;
        }

        locusResult[j] = locusFunction;
        total_cnts[j]  = overlapLen;
    }  // end Gene

    // Pick the most confidently
    if (result.size() > 1)
    {
        int           geneIdx    = 0;
        LocusFunction locus      = locusResult[geneIdx];
        int           overlapLen = total_cnts[geneIdx];
        bool          discard    = false;

        for (unsigned i = 1; i < result.size(); ++i)
        {
            LocusFunction tmpLocus = locusResult[i];
            int           tmpLen   = total_cnts[i];
            if (tmpLocus > locus)
            {
                // Find better gene
                geneIdx    = i;
                overlapLen = tmpLen;
                locus      = tmpLocus;
                discard    = false;
            }
            else if (tmpLocus == locus)
            {
                if (tmpLen > overlapLen)
                {
                    geneIdx    = i;
                    overlapLen = tmpLen;
                    discard    = false;
                }
                else if (tmpLen == overlapLen)
                {
                    discard = true;
                }
            }
        }

        if (discard)
            gidx = -1;
        else
            gidx = geneIdx;
    }

    // two gene has same overlap length, discard all of them
    if (gidx < 0)
    {
        overlap = 0;
        return { gidx, LocusFunction::INTERGENIC };
    }
    else
    {
        overlap = total_cnts[gidx];
        return { gidx, locusResult[gidx] };
    }
}

int TagReadsWithGeneExon::setAnnotation(BamRecord& record, const std::string& contig,
                                        string& gene, unsigned int* strand,
                                        unsigned int* locus)
{
    return annoByOverlapLen(record, contig, gene, strand, locus);
}

int TagReadsWithGeneExon::annoByOverlapLen(BamRecord& record, const std::string& contig,
                                           string& outGene, unsigned int* outStrand,
                                           unsigned int* outLocus)
{
    total_reads++;
    bool confidently = getQual(record) >= 255;
    if (confidently)
        map_reads++;

    std::vector< std::pair< int, int > > cigars = getCigar(record);

    // Change begin position from 0-based to 1-based
    int        beginPos   = getRefStart(record) + 1;
    int        queryBegin = beginPos;
    int        queryEnd   = queryBegin + getReferenceLength(cigars) - 1;
    MyInterval query_range{ queryBegin, queryEnd };
    std::vector< const GeneFromGTF* > result;

    // Skip genes with negative strand
    bool existsNegativeGene = false;
    bool readsNegative      = getNegativeStrand(record);
    if (mytrees.count(contig) != 0)
    {
        const auto& tree = mytrees.at(contig);
        const auto& res  = tree.query(query_range);
        for (const auto& o : res)
        {
            bool annoNegative = o.value.isNegativeStrand();
            if (readsNegative == annoNegative)
                result.push_back(&(o.value));
            else
                existsNegativeGene = true;
        }
    }

    if (existsNegativeGene)
        antisense_reads++;

    if (result.empty())
    {
        if (confidently)
        {
            intergenic_reads++;
        }
        *outLocus = 1;  // set intergenic locus
        return 0;
    }

    // Set annotations and record the metrics.
    std::vector< AlignmentBlock > alignmentBlock = getAlignmentBlocks(cigars, beginPos);
    int                           overlap;
    int                           idx;
    LocusFunction                 locus;
    if (multiMap)
        std::tie(idx, locus) = getLocusFunction(result, alignmentBlock, overlap);
    else
        std::tie(idx, locus) = getLocusFunction(result, alignmentBlock);

    *outLocus = int(locus);

    if (confidently)
    {
        if (locus == LocusFunction::CODING)
        {
            exonic_reads++;
            transcriptome_reads++;
        }
        else if (locus == LocusFunction::INTRONIC)
        {
            intronic_reads++;
            transcriptome_reads++;
        }
        else if (locus == LocusFunction::INTERGENIC)
            intergenic_reads++;
    }

    // Only dump gene name when locus is Exon or Intron
    if (locus != LocusFunction::INTERGENIC)
    {
        outGene    = result[idx]->getName();
        *outStrand = readsNegative
                         ? 2
                         : 1;  // same strand because we filter genes with negative strand
    }

    return 0;
}

int TagReadsWithGeneExon::setAnnotation(BamRecord& record, const std::string& contig,
                                        int& overlap, int& outLocus)
{
    std::vector< std::pair< int, int > > cigars = getCigar(record);

    // Change begin position from 0-based to 1-based
    int        beginPos   = getRefStart(record) + 1;
    int        queryBegin = beginPos;
    int        queryEnd   = queryBegin + getReferenceLength(cigars) - 1;
    MyInterval query_range{ queryBegin, queryEnd };
    std::vector< const GeneFromGTF* > result;

    // Skip genes with negative strand
    bool readsNegative = getNegativeStrand(record);
    if (mytrees.count(contig) != 0)
    {
        const auto& tree = mytrees.at(contig);
        const auto& res  = tree.query(query_range);
        for (const auto& o : res)
        {
            bool annoNegative = o.value.isNegativeStrand();
            if (readsNegative == annoNegative)
                result.push_back(&(o.value));
        }
    }

    if (result.empty())
    {
        outLocus = 1;  // set intergenic locus
        return 0;
    }

    // Set annotations and record the metrics.
    std::vector< AlignmentBlock > alignmentBlock = getAlignmentBlocks(cigars, beginPos);
    auto [idx, locus] = getLocusFunction(result, alignmentBlock, overlap);

    outLocus = int(locus);

    return 0;
}

int TagReadsWithGeneExon::makeOverlapDetector()
{
    libx::Timer                                                 load_timer;
    std::unordered_map< std::string, std::vector< GTFRecord > > gatherByGeneName;
    try
    {
        AnnoReader* annoReader = createAnnoReader(annotation_filename);
        annoReader->loadAnnoFile(gatherByGeneName);
        if (annoReader != nullptr)
        {
            delete annoReader;
            annoReader = nullptr;
        }
    }
    catch (AnnotationException& e)
    {
        spdlog::error(e.what());
        return -1;
    }
    catch (std::exception& e)
    {
        spdlog::error(e.what());
        return -1;
    }
    spdlog::info("loadGTFFile time(s):{}", load_timer.toc());

    int         numGene = 0;
    GeneBuilder geneBuilder;
    GTFIterator gtfIterator = gatherByGeneName.begin();

    geneName.push_back("");
    for (; gtfIterator != gatherByGeneName.end(); ++gtfIterator)
    {
        try
        {
            std::vector< GeneFromGTF > genes = geneBuilder.makeGene(gtfIterator);
            for (auto& gene : genes)
            {
                Node node;
                node.lower = gene.getStart();
                node.upper = gene.getEnd();
                node.value = gene;
                nodes.push_back(std::move(node));
                ++numGene;

                string gname = gene.getContig() + "_" + gene.getName();
                if (geneID.count(gname) == 0)
                {
                    geneID[gname] = geneName.size();
                    geneName.push_back(gene.getName());
                }
                geneEnd[gene.getContig()].push_back({ gname, gene.getEnd() });
            }
        }
        catch (AnnotationException& e)
        {
            spdlog::warn("{} {}", e.what(), "-- skipping");
        }
    }

    typedef pair< string, unsigned int > mypair;
    for (auto& [k, real] : geneEnd)
    {
        std::sort(real.begin(), real.end(),
                  [](mypair& t1, mypair& t2) { return t1.second < t2.second; });
    }

    for (auto& node : nodes)
        mytrees[node.value.getContig()].insert(node);
    spdlog::info("makeOverlapDetector time(s):{}", load_timer.toc());
    spdlog::info("Gene count:{}", numGene);
    if (numGene == 0)
    {
        spdlog::error("No valid gene found!");
        return -2;
    }

    return 0;
}

std::string TagReadsWithGeneExon::dumpMetrics()
{
    std::stringstream ss;
    ss << "## ANNOTATION METRICS\n";
    ss << "TOTAL_"
          "READS\tMAP\tEXONIC\tINTRONIC\tINTERGENIC\tTRANSCRIPTOME\tANTISENSE"
          "\n";
    ss << total_reads << "\t" << map_reads << "\t" << exonic_reads << "\t"
       << intronic_reads << "\t" << intergenic_reads << "\t" << transcriptome_reads
       << "\t" << antisense_reads << std::endl;
    ss.flags(std::ios::fixed);
    ss.precision(1);
    ss << 100.0 << "\t" << map_reads * 100.0 / total_reads << "\t"
       << exonic_reads * 100.0 / total_reads << "\t"
       << intronic_reads * 100.0 / total_reads << "\t"
       << intergenic_reads * 100.0 / total_reads << "\t"
       << transcriptome_reads * 100.0 / total_reads << "\t"
       << antisense_reads * 100.0 / total_reads << "\t" << std::endl;

    return ss.str();
}
