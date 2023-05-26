/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bamUtils.h"

using namespace std;

class Exon
{
public:
    Exon(int s, int e) : start(s), end(e) {}

    int start;
    int end;
};

class TranscriptFromGTF
{
public:
    TranscriptFromGTF(int transcriptionStart, int transcriptionEnd, int numExons,
                      std::string transcriptName, std::string transcriptID)
        : transcriptionStart(transcriptionStart), transcriptionEnd(transcriptionEnd),
          numExons(numExons), transcriptName(transcriptName), transcriptID(transcriptID)
    {
    }
    TranscriptFromGTF() {}

    void addExons(std::vector< Exon >& exons_);

    const std::string& getTranscriptName() const
    {
        return transcriptName;
    }
    int getTranscriptionStart() const
    {
        return transcriptionStart;
    }
    int getTranscriptionEnd() const
    {
        return transcriptionEnd;
    }
    const std::vector< Exon >& getExons() const
    {
        return exons;
    }

    void assignLocusFunction(int start, int len, pair< int, int >& max_cnts) const;
    bool inExon(int locus) const;

private:
    int         transcriptionStart;
    int         transcriptionEnd;
    int         numExons;
    std::string transcriptName;
    std::string transcriptID;

    std::vector< Exon > exons;
};

class GeneFromGTF
{
public:
    GeneFromGTF() {}
    GeneFromGTF(std::string& contig, std::string featureType, int start, int end,
                bool negativeStrand, std::string geneID, std::string& geneName)
        : contig(contig), featureType(featureType), start(start), end(end),
          negativeStrand(negativeStrand), geneID(geneID), geneName(geneName)
    {
    }

    bool isNegativeStrand() const
    {
        return negativeStrand;
    }
    const std::string& getContig() const
    {
        return contig;
    }
    const std::string& getName() const
    {
        return geneName;
    }
    int getStart() const
    {
        return start;
    }
    int getEnd() const
    {
        return end;
    }
    const std::string& getFeatureType() const
    {
        return featureType;
    }
    const std::unordered_map< std::string, TranscriptFromGTF >& getTranscripts() const
    {
        return transcripts;
    }

    TranscriptFromGTF* addTranscript(int transcriptionStart, int transcriptionEnd,
                                     int numExons, std::string transcriptName,
                                     std::string transcriptID);

private:
    std::string contig;
    std::string featureType;
    int         start;
    int         end;
    bool        negativeStrand;
    std::string geneID;
    std::string geneName;

    std::unordered_map< std::string, TranscriptFromGTF > transcripts;
};