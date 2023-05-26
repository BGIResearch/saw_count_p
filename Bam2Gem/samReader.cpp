/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
#include <filesystem>
namespace fs = std::filesystem;

#include <htslib/thread_pool.h>
#include <libx/System.hpp>
#include <libx/Timer.hpp>
#include <spdlog/spdlog.h>

#include "samReader.h"

// Class Implementations.
const int QUAL_THRESHOLD = 10;

static char SEP = '@';

// Common functions.
inline string parseQname(string& qname)
{
    string res;
    size_t pos0 = qname.find(SEP);
    res += qname.substr(0, pos0);
    size_t pos1 = qname.find(SEP, pos0 + 1);
    res += qname.substr(pos0 + 1, pos1 - pos0 - 1);
    return res;
}

// Class SamReader's functions.
SamReader::~SamReader()
{
    if (fp_)
    {
        Close();
    }
}

void check_index(const std::string& reads_path)
{
    // The htslib index data structure for our indexed BAM file. May be NULL if
    // no index was loaded.
    const int threads_num = 4;

    std::string bai_idx_path = reads_path + ".bai";
    std::string csi_idx_path = reads_path + ".csi";
    if (fs::exists(bai_idx_path) && libx::compareFileTime(reads_path, bai_idx_path))
        return;
    if (fs::exists(csi_idx_path) && libx::compareFileTime(reads_path, csi_idx_path))
        return;

    libx::Timer timer;
    // Try Create bai index first
    if (!fs::exists(bai_idx_path) || libx::compareFileTime(reads_path, bai_idx_path))
    {
        spdlog::debug("build index: {}", bai_idx_path);
        int ret =
            sam_index_build3(reads_path.c_str(), bai_idx_path.c_str(), 0, threads_num);
        spdlog::debug("build index time(s): {}", timer.toc());
        if (ret == 0)
            return;
    }

    // Bai-format limit < 512MB, try csi-format
    spdlog::debug("build index: {}", csi_idx_path);
    int ret = sam_index_build3(reads_path.c_str(), csi_idx_path.c_str(), 14, threads_num);
    spdlog::debug("build index time(s): {}", timer.toc());
    if (ret != 0)
        throw std::runtime_error(ERR_CODE + "003 " + "Failed to create index for "
                                 + reads_path);
}

std::unique_ptr< SamReader > SamReader::FromFile(const std::string& reads_path,
                                                 const int          block_size)
{
    htsFile* fp = hts_open(reads_path.c_str(), "r");
    if (!fp)
    {
        string errorMsg = "Could not open " + reads_path;
        spdlog::error(errorMsg);
        throw std::runtime_error(ERR_CODE + "002 " + errorMsg);
    }

    if (hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, block_size) != 0)
    {
        string errorMsg = "Failed to set HTS_OPT_BLOCK_SIZE as " + to_string(block_size);
        spdlog::error(errorMsg);
        throw std::runtime_error(ERR_CODE + "003 " + errorMsg);
    }

    bam_hdr_t* header = sam_hdr_read(fp);
    if (header == nullptr)
    {
        string errorMsg = "Failed to read bam header: " + reads_path;
        spdlog::error(errorMsg);
        throw std::runtime_error(ERR_CODE + "001 " + errorMsg);
    }

    // check_index(reads_path);
    hts_idx_t* idx = nullptr;
    // hts_idx_t* idx = sam_index_load(fp, fp->fn);
    // if (idx == nullptr)
    // {
    //     string errorMsg = "Failed to read bam index: " + reads_path;
    //     spdlog::error(errorMsg);
    //     throw std::runtime_error(ERR_CODE + "001 " + errorMsg);
    // }
    return std::unique_ptr< SamReader >(new SamReader(fp, header, idx));
}

int SamReader::QueryAll()
{
    bam1_t* b = bam_init1();
    // htsFile* out = hts_open("out.bam", "wb");
    // sam_hdr_write(out, header_);

    /*
    while (sam_read1(fp_, header_, b) >= 0)
    {
        if (b->core.qual < QUAL_THRESHOLD) continue;
        int r = sam_write1(out, header_, b);
    }
    */
    unordered_set< string > read_set;
    size_t                  total = 0, unique = 0;
    for (size_t i = 0; i < ref_.size(); ++i)
    {
        hts_itr_t* iter = sam_itr_queryi(idx_, i, 0, ref_[i].second);
        if (iter == nullptr)
        {
            cerr << "query unknown reference:" << ref_[i].first << endl;
            continue;
        }
        read_set.clear();
        while (sam_itr_next(fp_, iter, b) > 0)
        {
            // Filter reads which quality value less than quality threshold.
            if (b->core.qual < QUAL_THRESHOLD)
                continue;
            string qname  = bam_get_qname(b);
            string marker = parseQname(qname);
            if (read_set.count(marker) == 0)
            {
                ++unique;
                read_set.insert(marker);
                // sam_write1(out, header_, b);
            }
            ++total;
        }
    }
    cout << "Total Reads: " << total << "\nUnique Reads: " << unique << endl;
    cout << "Duplication rate: " << setiosflags(ios::fixed) << setprecision(2)
         << float(total - unique) * 100 / total << "%" << endl;

    bam_destroy1(b);
    // hts_close(out);
    // out = nullptr;

    return 0;
}

int SamReader::QueryOne(BamRecord b)
{
    for (size_t i = 0; i < ref_.size(); ++i)
    {
        hts_itr_t* iter = sam_itr_queryi(idx_, i, 0, ref_[i].second);
        if (iter == nullptr)
        {
            spdlog::warn("query unknown reference:{}", ref_[i].first);
            continue;
        }
        if (sam_itr_next(fp_, iter, b) > 0)
        {
            break;
        }
    }
    return 0;
}

int SamReader::QueryAll(BamRecord b)
{
    if (sam_read1(fp_, header_, b) >= 0)
    {
        return true;
    }
    return false;
}

std::vector< std::pair< std::string, uint32 > > SamReader::getContigs()
{
    return ref_;
}

bool SamReader::QueryByContig(int tid)
{
    iter_ = sam_itr_queryi(idx_, tid, 0, ref_[tid].second);
    if (iter_ == nullptr || iter_->finished)
    {
        // spdlog::debug("No reads for ref:{}", ref_[tid].first);
        return false;
    }
    return true;
}

bool SamReader::QueryByContigBE(int tid, const int beg, const int end, hts_itr_t*& iter)
{
    iter = sam_itr_queryi(idx_, tid, beg, end);
    if (iter == nullptr)
    {
        cerr << "query unknown reference:" << ref_[tid].first << endl;
        return false;
    }

    return true;
}

bool SamReader::next(BamRecord b)
{
    if (iter_ != nullptr && sam_itr_next(fp_, iter_, b) > 0)
    {
        return true;
    }
    return false;
}

std::string SamReader::refName(BamRecord b)
{
    const int tid = b->core.tid;
    if (tid >= 0)
        return header_->target_name[tid];
    else
        return "";
}

bool SamReader::next(BamRecord b, hts_itr_t*& iter)
{
    int ret;
    if ((ret = sam_itr_next(fp_, iter, b)) > 0)
    {
        return true;
    }
    return false;
}

SamReader::SamReader(htsFile* fp, bam_hdr_t* header, hts_idx_t* idx)
    : fp_(fp), header_(header), idx_(idx)
{
    iter_ = nullptr;
    for (int i = 0; i < header->n_targets; ++i)
        ref_.push_back({ header_->target_name[i], header_->target_len[i] });
}

int SamReader::Close()
{
    if (fp_ != nullptr)
    {
        hts_close(fp_);
        fp_ = nullptr;
    }
    if (idx_ != nullptr)
    {
        hts_idx_destroy(idx_);
        idx_ = nullptr;
    }
    if (iter_ != nullptr)
    {
        hts_itr_destroy(iter_);
        iter_ = nullptr;
    }
    if (header_ != nullptr)
    {
        bam_hdr_destroy(header_);
        header_ = nullptr;
    }
    return 0;
}

BamHeader SamReader::getHeader()
{
    return header_;
}

int SamReader::setThreadPool(htsThreadPool* p)
{
    return hts_set_thread_pool(fp_, p);
}

int SamReader::setThreads(int n)
{
    return hts_set_threads(fp_, n);
}

std::map< std::string, uint64_t > SamReader::getReadsByContig()
{
    std::map< std::string, uint64_t > res;

    for (int i = 0; i < sam_hdr_nref(header_); ++i)
    {
        string   contig = sam_hdr_tid2name(header_, i);
        uint64_t u, v;
        hts_idx_get_stat(idx_, i, &u, &v);
        res[contig] = u + v;
    }

    return res;
}

MergeSamReader::MergeSamReader(std::vector< std::string >& filenames)
{
    for (auto& filename : filenames)
    {
        auto reader = SamReader::FromFile(filename);
        readers.emplace_back(std::move(reader));
    }

    finished       = false;
    readsExhausted = false;
    beginRead      = false;
    // TODO: make sure all the bam files have same contigs
}

MergeSamReader::~MergeSamReader()
{
    for (auto& reader : readers)
    {
        reader->Close();
    }
}

bool MergeSamReader::QueryByContig(int tid)
{
    bool res = false;
    for (unsigned int i = 0; i < readers.size(); ++i)
    {
        if (readers[i]->QueryByContig(tid))
            res = true;
    }
    // reset status of reads exhausted
    readsExhausted = false;

    // begin read bam files
    beginRead = true;

    return res;
}

bool MergeSamReader::next(BamData& b)
{
    while (true)
    {
        if (bamQueue.try_dequeue(b))
        {
            return true;
        }
        if (readsExhausted)
            break;
    }
    return false;
}

void MergeSamReader::daemon(unsigned int threadNum)
{
    libx::Timer       timer("s");
    vector< BamData > records(readers.size());
    vector< bool >    recordStatus(readers.size(), false);

    htsThreadPool p = { NULL, 0 };
    p.pool          = hts_tpool_init(threadNum);
    if (!p.pool)
    {
        spdlog::warn("fail init hts thread pool for reading");
    }
    else
    {
        for (auto& reader : readers)
            reader->setThreadPool(&p);
    }

    auto         contigs = readers.front()->getContigs();
    unsigned int chrID   = 0;

    int currStart = INT_MAX;
    int currPos   = -1;

    bool chrExhausted = false;

    for (unsigned int i = 0; i < readers.size(); ++i)
    {
        readers[i]->next(records[i]._record);
    }
    unsigned long reads = 0;
    while (!finished)
    {
        // incase use too much memory
        // if (bamQueue.size_approx() > 1000)
        // {
        //     std::this_thread::yield();
        //     continue;
        // }

        // sort bam by mapping position
        // QueryByContig() and next() must be mutually exclusive
        currStart = INT_MAX;
        currPos   = -1;
        for (unsigned int i = 0; i < readers.size(); ++i)
        {
            if (!recordStatus[i])
                continue;
            int refStart = getRefStart(records[i]._record);
            if (refStart <= currStart)
            {
                currStart = refStart;
                currPos   = i;
            }
        }

        if (currPos == -1)
        {
            spdlog::debug("processed reads: {} time(s): {}", reads, timer.toc());
            if (chrID == contigs.size())
                break;
            while (chrID < contigs.size())
            {
                if (!QueryByContig(chrID))
                {
                    chrID++;
                    if (chrID == contigs.size())
                    {
                        chrExhausted = true;
                        break;
                    }
                }
                else
                {
                    chrID++;
                    break;
                }
            }
            if (chrExhausted)
                break;

            for (unsigned int i = 0; i < readers.size(); ++i)
            {
                if (readers[i]->next(records[i]._record))
                    recordStatus[i] = true;
                else
                    recordStatus[i] = false;
            }

            continue;
        }

        bamQueue.wait_enqueue(records[currPos]);
        reads++;

        if (!readers[currPos]->next(records[currPos]._record))
        {
            // records[currPos] = nullptr;
            recordStatus[currPos] = false;
        }
    }
    spdlog::debug("reader finished");

    // reads has been exhausted
    readsExhausted = true;

    if (p.pool)
    {
        hts_tpool_destroy(p.pool);
    }

    // for (auto& record : records)
    //     destroyBamRecord(record);
    // destroyBamRecord(currRecord);
}

void MergeSamReader::setFinish()
{
    finished = true;
}