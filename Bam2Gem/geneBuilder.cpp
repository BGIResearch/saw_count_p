/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include <climits>
#include <iostream>

#include <algorithm>
#include <set>

#include "annotationException.h"
#include "geneBuilder.h"

// For logging format
string geneIDsToString(set< string >& ss, string sep = ", ")
{
    string res;
    for (auto& s : ss)
    {
        if (!res.empty())
            res += sep;
        res += s;
    }
    return res;
}

vector< GeneFromGTF > GeneBuilder::makeGene(GTFIterator& iter)
{
    vector< GTFRecord > gtfRecords = iter->second;
    // Check is there same gene in different chromosomes
    if (geneInDiffChrs(gtfRecords))
    {
        vector< GeneFromGTF > res;
        // Grouping input gtf records into multi fragments
        vector< GTFRecord > temp;
        for (auto& record : gtfRecords)
        {
            if (!temp.empty() && record.getContig() != temp.front().getContig())
            {
                res.push_back(makeGeneWithTranscripts(temp));
                temp.clear();
            }
            temp.push_back(record);
        }
        if (!temp.empty())
        {
            res.push_back(makeGeneWithTranscripts(temp));
            temp.clear();
        }
        return res;
    }
    else
    {
        return { makeGeneWithTranscripts(gtfRecords) };
    }
}

bool GeneBuilder::geneInDiffChrs(vector< GTFRecord >& gtfRecords)
{
    string contig(gtfRecords[0].getContig());
    for (auto& record : gtfRecords)
    {
        if (record.getContig() != contig)
            return true;
    }
    return false;
}

GeneFromGTF GeneBuilder::makeGeneWithTranscripts(vector< GTFRecord >& gtfRecords)
{
    GeneFromGTF gene = makeGene(gtfRecords);

    GTFMapForTranscript gftRecordsByTranscript = gatherByTransciptID(gtfRecords);
    for (auto& [_, gftRecords] : gftRecordsByTranscript)
        addTranscriptToGene(gene, gftRecords);

    if (gene.getTranscripts().empty())
    {
        throw AnnotationException("No transcript for gene " + gene.getName());
    }
    return gene;
}

GTFMapForTranscript GeneBuilder::gatherByTransciptID(vector< GTFRecord >& gtfRecords)
{
    GTFMapForTranscript res;
    for (auto& record : gtfRecords)
    {
        // Skip records with feature type "gene" when makeing transcripts.
        if (record.getFeatureType() == "gene")
            continue;
        if (record.getTranscriptID() == "")
        {
            throw AnnotationException("Record does not have transcriptID for gene "
                                      + record.getGeneName());
        }
        res[record.getTranscriptID()].push_back(record);
    }
    return res;
}

GeneFromGTF GeneBuilder::makeGene(vector< GTFRecord >& gtfRecords)
{
    auto& geneRecord = gtfRecords[0];

    // Figure out the extend of the gene.
    int           start = INT_MAX;
    int           end   = INT_MIN;
    set< string > geneIds;
    uint32_t      geneNum = 0;
    for (auto& record : gtfRecords)
    {
        start = min(start, record.getStart());
        end   = max(end, record.getEnd());
        geneIds.insert(record.getGeneID());
        if (record.getFeatureType() == "gene")
        {
            geneNum++;
        }
    }

    if (geneIds.size() != geneNum)
    {
        string error = "Multiple gene IDs for gene " + geneRecord.getGeneName() + ": "
                       + geneIDsToString(geneIds);
        throw AnnotationException(error);
    }

    GeneFromGTF gene(geneRecord.getContig(), geneRecord.getFeatureType(), start, end,
                     geneRecord.isNegativeStrand(), geneRecord.getGeneID(),
                     geneRecord.getGeneName());

    for (auto& record : gtfRecords)
        validateGTFRecord(record, gene);

    return gene;
}

void GeneBuilder::validateGTFRecord(GTFRecord& record, GeneFromGTF& gene)
{
    if (gene.isNegativeStrand() != record.isNegativeStrand())
    {
        throw AnnotationException("Strand disagreement for gene " + gene.getName());
    }
}

void GeneBuilder::addTranscriptToGene(GeneFromGTF& gene, vector< GTFRecord >& gtfRecords)
{
    auto&  geneName              = gene.getName();
    auto&  lineOne               = gtfRecords[0];
    auto&  transcriptName        = lineOne.getTranscriptName();
    auto&  transcriptID          = lineOne.getTranscriptID();
    string transcriptDescription = geneName + ":" + transcriptName;

    vector< Exon > exons;
    int            transcriptionStart = INT_MAX;
    int            transcriptionEnd   = INT_MIN;
    for (auto& record : gtfRecords)
    {
        string featureType = record.getFeatureType();
        int    start       = record.getStart();
        int    end         = record.getEnd();
        if (featureType == "exon")
        {
            Exon e(start, end);
            exons.push_back(std::move(e));
            transcriptionStart = min(transcriptionStart, start);
            transcriptionEnd   = max(transcriptionEnd, end);
        }
    }

    // Sort exons by start
    sort(exons.begin(), exons.end(),
         [](Exon& e1, Exon& e2) { return (e1.start < e2.start); });
    // Validate exons
    for (size_t i = 0; i < exons.size(); ++i)
    {
        auto& e = exons[i];
        if (e.start > e.end)
        {
            throw AnnotationException("Exon has 0 or negative extent for "
                                      + transcriptDescription);
        }
        if (i > 0 && exons[i - 1].end > e.start)
        {
            throw AnnotationException("Exons overlap for " + transcriptDescription);
        }
    }

    TranscriptFromGTF* transcript = gene.addTranscript(
        transcriptionStart, transcriptionEnd, exons.size(), transcriptName, transcriptID);
    transcript->addExons(exons);
}
