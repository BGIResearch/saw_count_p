/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <string>
using std::string;

class GTFRecord
{
public:
    GTFRecord(string& contig, string& featureType, int start, int end,
              bool negativeStrand, string& geneID, string& geneName, string& transcriptID,
              string& transcriptName)
        : contig(contig), featureType(featureType), start(start), end(end),
          negativeStrand(negativeStrand), geneID(geneID), geneName(geneName),
          transcriptID(transcriptID), transcriptName(transcriptName)
    {
    }

    // Get methods
    string& getGeneID()
    {
        return geneID;
    }
    string& getGeneName()
    {
        return geneName;
    }
    string& getTranscriptID()
    {
        return transcriptID;
    }
    string& getTranscriptName()
    {
        return transcriptName;
    }
    string& getFeatureType()
    {
        return featureType;
    }
    int getStart()
    {
        return start;
    }
    int getEnd()
    {
        return end;
    }
    string& getContig()
    {
        return contig;
    }
    bool isNegativeStrand()
    {
        return negativeStrand;
    }

private:
    string contig;
    string featureType;
    int    start, end;
    bool   negativeStrand;

    string geneID;
    string geneName;

    string transcriptID;
    string transcriptName;
};
