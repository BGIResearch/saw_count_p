/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include <iostream>

#include <spdlog/spdlog.h>

#include "annotationException.h"
#include "geneFromGTF.h"

inline bool inRange(int start, int end, int locus)
{
    return start <= locus && locus <= end;
}

TranscriptFromGTF* GeneFromGTF::addTranscript(int transcriptionStart,
                                              int transcriptionEnd, int numExons,
                                              std::string transcriptName,
                                              std::string transcriptID)
{
    // Keep transcripts with same name
    if (transcripts.count(transcriptID) != 0)
    {
        std::string error =
            "Transcript appears more than once for " + geneName + " : " + transcriptID;
        throw AnnotationException(error);
    }
    else
    {
        TranscriptFromGTF transcript(transcriptionStart, transcriptionEnd, numExons,
                                     transcriptName, transcriptID);
        transcripts[transcriptID] = transcript;
        return &(transcripts[transcriptID]);
    }
}

void TranscriptFromGTF::addExons(std::vector< Exon >& exons_)
{
    exons.assign(exons_.begin(), exons_.end());
}

void TranscriptFromGTF::assignLocusFunction(int start, int len,
                                            pair< int, int >& max_cnts) const
{
    int begin = std::max(start, transcriptionStart);
    int end   = std::min(int(start + len - 1), transcriptionEnd);
    if (begin > end)
        return;
    int exon_cnts = 0, intro_cnts = 0;
    for (int i = begin; i <= end; ++i)
    {
        if (inExon(i))
            ++exon_cnts;
        else
            ++intro_cnts;
    }
    if (exon_cnts > max_cnts.first)
        max_cnts = { exon_cnts, intro_cnts };
    else if (exon_cnts == max_cnts.first && intro_cnts > max_cnts.second)
        max_cnts.second = intro_cnts;
}

bool TranscriptFromGTF::inExon(int locus) const
{
    for (auto& exon : exons)
    {
        if (exon.start > locus)
            return false;
        if (inRange(exon.start, exon.end, locus))
            return true;
    }
    return false;
}
