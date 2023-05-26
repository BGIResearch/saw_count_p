/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

typedef bam1_t*    BamRecord;
typedef bam_hdr_t* BamHeader;

// 0->None
// 1->INTERGENIC->2
// 2->INTRONIC->1
// 3->EXONIC->0
static const int LocusTransfer[] = { 0, 2, 1, 0 };

/**
 * Represents the contiguous alignment of a subset of read bases to a reference
 * sequence. Simply put an alignment block tells you that read bases from
 * readStart are aligned to the reference (matching or mismatching) from
 * referenceStart for length bases.
 */
class AlignmentBlock
{
public:
    /** Constructs a new alignment block with the supplied read and ref starts
     * and length. */
    AlignmentBlock() {}

    AlignmentBlock(int readStart_, int referenceStart_, int length_)
        : readStart(readStart_), referenceStart(referenceStart_), length(length_)
    {
    }

    /** The first, 1-based, base in the read that is aligned to the reference
     * reference. */
    int getReadStart()
    {
        return readStart;
    }

    /** The first, 1-based, position in the reference to which the read is
     * aligned. */
    int getReferenceStart()
    {
        return referenceStart;
    }

    /** The number of contiguous bases aligned to the reference. */
    int getLength()
    {
        return length;
    }

private:
    int readStart;
    int referenceStart;
    int length;
};

/**
 * Given a Cigar, Returns blocks of the sequence that have been aligned directly
 * to the reference sequence. Note that clipped portions, and inserted and
 * deleted bases (vs. the reference) are not represented in the alignment
 * blocks.
 *
 * @param cigar          The cigar containing the alignment information
 * @param alignmentStart The start (1-based) of the alignment
 * @param cigarTypeName  The type of cigar passed - for error logging.
 * @return List of alignment blocks
 */
std::vector< AlignmentBlock >
getAlignmentBlocks(std::vector< std::pair< int, int > >& cigars, int alignmentStart);

int getReferenceLength(std::vector< std::pair< int, int > >& cigars);

int getQName(BamRecord b, std::string& qname);

int getQual(BamRecord b);

std::vector< std::pair< int, int > > getCigar(BamRecord b);

int getRefStart(BamRecord b);

bool getNegativeStrand(BamRecord b);

int updateStrTags(BamRecord& b, std::string tag, std::string data);

int updateIntTags(BamRecord& b, std::string tag, unsigned int data);

BamRecord createBamRecord();

int destroyBamRecord(BamRecord b);

bool compareBamRecord(BamRecord b1, BamRecord b2);

bool getTag(BamRecord b, const char tag[2], std::string& value);

bool getTagInt(BamRecord b, const char tag[2], int& value);

void setQcFail(BamRecord b);

void setDuplication(BamRecord b);

bool getQcFail(BamRecord b);

bool getDuplication(BamRecord b);

int getMarker(BamRecord b, std::string& marker);

class BamData
{
public:
    BamData();
    ~BamData();

    BamData& operator=(BamData& other);

public:
    BamRecord _record;
};