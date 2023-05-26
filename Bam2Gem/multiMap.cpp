/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include <regex>

#include <spdlog/spdlog.h>

#include "multiMap.h"

#include "libx/System.hpp"
#include "libx/ThreadPool.hpp"
#include "libx/Timer.hpp"

MultiMap::MultiMap(vector< string >& bamFiles, uint32 threads)
    : bamFiles(bamFiles), threads(threads), qnameStrMode(false)
{
    int len      = bamFiles.size();
    blockSize    = std::max(_DEFAULT_HTS_BLOCK_SIZE / len, 1024 * 1024);
    threadPerBam = std::max(int(threads / len), 1);

    strUmaps.resize(len);
    intUmaps.resize(len);
}

void MultiMap::collectMultiMap(int idx, TagReadsWithGeneExon* anno)
{
    libx::Timer timer;

    auto reader = SamReader::FromFile(bamFiles[idx], blockSize);
    auto contigs = reader->getContigs();
    reader->setThreads(threadPerBam);
    BamRecord read    = createBamRecord();

    auto& smap = strUmaps[idx];
    auto& imap = intUmaps[idx];

    int    nh, hi;
    string qname;
    
    while (reader->QueryAll(read))
    {
        string ctg = contigs[read->core.tid].first;
        getTagInt(read, "NH", nh);
        if (nh <= 1)
            continue;

        int overlap = 0, locus = 0;
        anno->setAnnotation(read, ctg, overlap, locus);

        getQName(read, qname);
        getTagInt(read, "HI", hi);
        if (qnameStrMode)
        {
            if (smap.count(qname) == 0)
            {
                ReadInfo mm(overlap, locus, nh, hi);
                smap[qname] = mm;
            }
            else
            {
                smap[qname].update(overlap, locus, hi);
            }
        }
        else
        {
            uint64 key = encodeQnameV2(qname);
            if (imap.count(key) == 0)
            {
                ReadInfo mm(overlap, locus, nh, hi);
                imap[key] = mm;
            }
            else
            {
                imap[key].update(overlap, locus, hi);
            }
        }
    }
       
    destroyBamRecord(read);

    uint64 mmReads = qnameStrMode ? smap.size() : imap.size();
    spdlog::debug("finish create map for bam: {} mm reads: {} time(s): {}", idx, mmReads,
                  timer.toc());
}

bool MultiMap::createMap(TagReadsWithGeneExon* anno)
{
    libx::Timer timer;
    bool        res = true;
    spdlog::debug("memory(GB) before create map of multi-map: {:.3f}",
                  libx::getSelfMemory());

    qnameStrMode = isStrQname();
    // qnameStrMode = true;
    spdlog::debug("use str qname: {}", qnameStrMode);

    if (!tryCreateMap(anno))
    {
        if (!qnameStrMode)
        {
            spdlog::debug("switch to str qname");
            intUmaps.clear();
            qnameStrMode = true;
            if (!tryCreateMap(anno))
                res = false;
        }
        else
            res = false;
    }

    spdlog::debug("memory(GB) after create map of multi-map: {:.3f} time(s): {}",
                  libx::getSelfMemory(), timer.toc());

    return res;
}

bool MultiMap::tryCreateMap(TagReadsWithGeneExon* anno)
{
    try
    {
        libx::ThreadPool         thpool(threads);
        vector< future< void > > results;
        for (uint32 idx = 0; idx < bamFiles.size(); ++idx)
        {
            results.emplace_back(thpool.commit(
                std::bind(&MultiMap::collectMultiMap, this, idx, std::ref(anno))));
        }
        for (auto& result : results)
            result.get();
    }
    catch (std::exception& e)
    {
        cerr << "error when try to create map for multi-map reads: " << e.what() << endl;
        return false;
    }

    return true;
}

int MultiMap::search(int idx, string& qname)
{
    int res = 0;
    if (qnameStrMode)
    {
        if (strUmaps[idx].count(qname))
            res = strUmaps[idx][qname].hi;
    }
    else
    {
        uint64 key = encodeQnameV2(qname);
        if (intUmaps[idx].count(key))
            res = intUmaps[idx][key].hi;
    }
    return res;
}

bool MultiMap::isStrQname()
{
    bool      res       = false;
    BamRecord bamRecord = createBamRecord();
    for (auto& bam : bamFiles)
    {
        auto reader = SamReader::FromFile(bam, blockSize);
        reader->QueryAll(bamRecord);
        string qname;
        getQName(bamRecord, qname);
        if (encodeQname(qname) == 0)
            res = true;
    }
    destroyBamRecord(bamRecord);

    return res;
}

uint64 MultiMap::encodeQname(string& qname)
{
    uint64 res = 0;
    regex  reg(R"(^(\w+)L(\d{1})C(\d{3})R(\d{3})(\d{8})$)");
    smatch match;

    if (regex_match(qname, match, reg))
    {
        string sn(match[1]);
        string snStr;
        if (snMap.count(sn))
        {
            snStr = snMap[sn];
        }
        else
        {
            snStr     = std::to_string(snMap.size());
            snMap[sn] = snStr;
        }
        string resStr = snStr;
        for (int i = 2; i <= 6; ++i)
            resStr += match[i];

        res = std::stoul(resStr);
    }

    return res;
}

uint64 MultiMap::encodeQnameV2(string& qname)
{
    uint64 res = 0;

    int    len = qname.size();
    string dnb = qname.substr(len - 8);
    string fovc, fovr, lane, sn;
    fovr = qname.substr(len - 8 - 3, 3);
    fovc = qname.substr(len - 8 - 4 - 3, 3);
    lane = qname.substr(len - 8 - 4 - 4 - 1, 1);
    sn   = qname.substr(0, len - 8 - 4 - 4 - 2);

    string snStr;
    if (snMap.count(sn))
    {
        snStr = snMap[sn];
    }
    else
    {
        snStr     = std::to_string(snMap.size());
        snMap[sn] = snStr;
    }

    string resStr = snStr + lane + fovc + fovr + dnb;
    res           = std::stoul(resStr);

    return res;
}