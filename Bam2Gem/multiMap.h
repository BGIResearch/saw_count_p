/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

#include "bamUtils.h"
#include "samReader.h"
#include "tagReadsWithGeneExon.h"
#include "utils.h"

// read info for multi-mapping correction
struct ReadInfo
{
    ReadInfo() {}
    ReadInfo(int overlap, int locus, int nh, int hi)
        : overlap(overlap), locus(locus), nh(nh), hi(hi)
    {
        // locus should be 1/2/3
        // represent intergenic/intronic/exonic
        // intergenic should be discard
        if (locus == 1)
            hi = 0;
    }
    void update(int _overlap, int _locus, int _hi)
    {
        // Got a lower level locus or it is intergenic, do nothing
        if (_locus < locus || _locus == 1)
            return;
        // Got a higher level locus, replace old data
        if (_locus > locus)
        {
            overlap = _overlap;
            locus   = _locus;
            hi      = _hi;
        }
        else
        {
            // Compare by length, and pick the read with longest overlap
            if (_overlap > overlap)
            {
                overlap = _overlap;
                locus   = _locus;
                hi      = _hi;
            }
            else if (_overlap == overlap)
            {
                hi = 0;
            }
        }
    }

    uint8 overlap;
    uint8 locus;
    uint8 nh;
    uint8 hi;
};

typedef unordered_map< string, ReadInfo > strMap;
typedef unordered_map< uint64, ReadInfo > intMap;

class MultiMap
{
public:
    MultiMap(vector< string >& bamFiles, uint32 threads);
    ~MultiMap() {}

public:
    // travel bam files and create map by qname
    // return false if failed
    bool createMap(TagReadsWithGeneExon* tagReadsWithGeneExon);
    // search hi index by qname and position
    // return 0 if qname not exists
    int search(int idx, string& qname);

private:
    bool tryCreateMap(TagReadsWithGeneExon* tagReadsWithGeneExon);
    void collectMultiMap(int idx, TagReadsWithGeneExon* anno);

    // return true for using str mode
    // return false for using uint64 mode
    bool isStrQname();

    // encode qname from string to uint64
    // return 0 if qname is invalid
    uint64 encodeQname(string& qname);
    // speed up and stricter format requirements
    uint64 encodeQnameV2(string& qname);

private:
    // parameters
    vector< string > bamFiles;
    uint32           threads;

    bool                            qnameStrMode;
    unordered_map< string, string > snMap;  // encode slide number in qname prefix
    vector< strMap >                strUmaps;
    vector< intMap >                intUmaps;

    int blockSize;
    int threadPerBam;
};