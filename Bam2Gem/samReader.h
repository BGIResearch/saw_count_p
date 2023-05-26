/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <readerwritercircularbuffer.h>
using namespace moodycamel;
#include <libx/Timer.hpp>

#include "bamUtils.h"
#include "utils.h"

class SamReader
{
public:
    // Creates a new SamReader reading from the BAM file reads_path
    static std::unique_ptr< SamReader >
    FromFile(const std::string& reads_path,
             const int          block_size = _DEFAULT_HTS_BLOCK_SIZE);

    ~SamReader();

    // DIsable assignment/copy operations.
    SamReader(const SamReader& other) = delete;
    SamReader& operator=(const SamReader&) = delete;

    // Gets all of the reads that overlap any bases in range.
    int QueryAll();

    // Query reads treat the bam as a whole contig
    int QueryAll(BamRecord b);

    // Get first bam record in order
    int QueryOne(BamRecord b);

    // Close the underlying resource descriptors.
    int Close();

    // Return <chromosome,index> pairs.
    std::vector< std::pair< std::string, uint32 > > getContigs();

    // Return header for construct other instance.
    BamHeader getHeader();

    // Query records by a given contig.
    bool QueryByContig(int tid);

    bool QueryByContigBE(int tid, const int beg, const int end, hts_itr_t*& iter);

    // Iterate all records from query range.
    // NOTICE: must call next() after call QueryByContig()
    bool next(BamRecord b);

    bool next(BamRecord b, hts_itr_t*& iter);

    int setThreadPool(htsThreadPool* p);
    int setThreads(int n);

    // Parse the contig name of read
    std::string refName(BamRecord b);

    // Parse real reads number from header and index file
    std::map< std::string, uint64_t > getReadsByContig();

private:
    // Private constructor; use FromFile to safely create a SamReader from a
    // file.
    SamReader(htsFile* fp, bam_hdr_t* header, hts_idx_t* idx);

    // A pointer to the htslib file used to access the SAM/BAM data.
    htsFile* fp_;

    // A htslib header data structure obtained by parsing the header of this
    // BAM.
    bam_hdr_t* header_;

    // The htslib index data structure for our indexed BAM file.
    // May be NULL if no index was loaded.
    hts_idx_t* idx_;

    // Store reference name and length from header.
    std::vector< std::pair< std::string, uint32 > > ref_;

    // A pointer to the query result.
    hts_itr_t* iter_;

    double nextTimes;
};

class MergeSamReader
{
public:
    MergeSamReader(std::vector< std::string >& filenames);
    ~MergeSamReader();

public:
    bool QueryByContig(int tid);
    bool next(BamData& b);
    void setFinish();
    void daemon(unsigned int threadNum);

private:
private:
    std::vector< std::unique_ptr< SamReader > > readers;

    BlockingReaderWriterCircularBuffer< BamData > bamQueue{ 1000 };
    std::atomic< bool >                           finished;
    std::atomic< bool >                           readsExhausted;
    std::atomic< bool >                           beginRead;
};