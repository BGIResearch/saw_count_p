/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include "bamUtils.h"

/*
cigars's op:
BAM_CMATCH      0	M
BAM_CINS        1	I
BAM_CDEL        2	D
BAM_CREF_SKIP   3	N
BAM_CSOFT_CLIP  4	S
BAM_CHARD_CLIP  5	H
BAM_CPAD        6	P
BAM_CEQUAL      7	EQ/=
BAM_CDIFF       8	X
BAM_CBACK       9
*/
std::vector< AlignmentBlock >
getAlignmentBlocks(std::vector< std::pair< int, int > >& cigars, int alignmentStart)
{
    if (cigars.empty())
        return {};
    std::vector< AlignmentBlock > alignmentBlocks;
    int                           readBase = 1;
    int                           refBase  = alignmentStart;
    int                           len;
    for (auto& c : cigars)
    {
        switch (c.first)
        {
        case BAM_CHARD_CLIP:
            break;  // ignore hard clips
        case BAM_CPAD:
            break;  // ignore pads
        case BAM_CSOFT_CLIP:
            readBase += c.second;
            break;  // soft clip read bases
        case BAM_CREF_SKIP:
            refBase += c.second;
            break;  // reference skip
        case BAM_CDEL:
            refBase += c.second;
            break;
        case BAM_CINS:
            readBase += c.second;
            break;
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            len = c.second;
            alignmentBlocks.push_back(AlignmentBlock(readBase, refBase, len));
            readBase += len;
            refBase += len;
            break;
        default:
            break;
            // throw std::exception("Case statement didn't deal op: " +
            // c.operation);
        }
    }
    return alignmentBlocks;
}

int getReferenceLength(std::vector< std::pair< int, int > >& cigars)
{
    int length = 0;
    for (auto& c : cigars)
    {
        switch (c.first)
        {
        case BAM_CMATCH:
        case BAM_CDEL:
        case BAM_CREF_SKIP:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            length += c.second;
            break;
        default:
            break;
        }
    }
    return length;
}

int getQName(BamRecord b, std::string& qname)
{
    qname = bam_get_qname(b);
    return 0;
}

int getQual(BamRecord b)
{
    return b->core.qual;
}

std::vector< std::pair< int, int > > getCigar(BamRecord b)
{
    std::vector< std::pair< int, int > > cigars(b->core.n_cigar);

    uint32_t* cigar = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i)
    {
        int op    = bam_cigar_op(cigar[i]);
        int len   = bam_cigar_oplen(cigar[i]);
        cigars[i] = { op, len };
    }
    return cigars;
}

int getRefStart(BamRecord b)
{
    return b->core.pos;
}

bool getNegativeStrand(BamRecord b)
{
    return bam_is_rev(b);
}

int updateStrTags(BamRecord& b, std::string tag, std::string data)
{
    return bam_aux_update_str(b, tag.c_str(), data.size() + 1, data.c_str());
}

int updateIntTags(BamRecord& b, std::string tag, unsigned int data)
{
    return bam_aux_update_int(b, tag.c_str(), data);
}

BamRecord createBamRecord()
{
    return bam_init1();
}

int destroyBamRecord(BamRecord b)
{
    bam_destroy1(b);
    return 0;
}

bool compareBamRecord(BamRecord b1, BamRecord b2)
{
    // if (getContig(b1) != getContig(b2))
    //     return false;
    if (getRefStart(b1) != getRefStart(b2))
        return false;
    if (getNegativeStrand(b1) != getNegativeStrand(b2))
        return false;
    if (b1->core.n_cigar != b2->core.n_cigar)
        return false;
    else
    {
        uint32_t* cigar1 = bam_get_cigar(b1);
        uint32_t* cigar2 = bam_get_cigar(b2);
        for (uint32_t i = 0; i < b1->core.n_cigar; ++i)
        {
            int op1  = bam_cigar_op(cigar1[i]);
            int len1 = bam_cigar_oplen(cigar1[i]);
            int op2  = bam_cigar_op(cigar2[i]);
            int len2 = bam_cigar_oplen(cigar2[i]);
            if (op1 != op2 || len1 != len2)
                return false;
        }
    }
    return true;
}

// Used for deduplication when no umi found
int getMarker(BamRecord b, std::string& marker)
{
    if (b->core.isize < 0)
    {
        marker =
            std::to_string(b->core.mpos - b->core.isize) + std::to_string(b->core.mpos);
    }
    else
    {
        marker =
            std::to_string(b->core.pos) + std::to_string(b->core.pos + b->core.isize);
    }
    return 0;
}

bool getTag(BamRecord b, const char tag[2], std::string& value)
{
    uint8_t* data = bam_aux_get(b, tag);
    if (data == NULL)
        return false;
    value = bam_aux2Z(data);
    return true;
}

bool getTagInt(BamRecord b, const char tag[2], int& value)
{
    uint8_t* data = bam_aux_get(b, tag);
    if (data == NULL)
        return false;
    value = bam_aux2i(data);
    return true;
}

void setQcFail(BamRecord b)
{
    b->core.flag |= BAM_FQCFAIL;
}

void setDuplication(BamRecord b)
{
    b->core.flag |= BAM_FDUP;
}

bool getQcFail(BamRecord b)
{
    return b->core.flag & BAM_FQCFAIL;
}

bool getDuplication(BamRecord b)
{
    return b->core.flag & BAM_FDUP;
}

BamData::BamData()
{
    _record = createBamRecord();
}

BamData::~BamData()
{
    destroyBamRecord(_record);
}

// void BamData::copy(BamRecord record)
// {
//     bam_copy1(_record, record);
// }

BamData& BamData::operator=(BamData& other)
{
    [[maybe_unused]] const auto _ = bam_copy1(_record, other._record);
    return *this;
}