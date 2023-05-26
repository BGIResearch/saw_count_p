/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#pragma once

#include <filesystem>  // C++17 only
#include <sstream>
#include <tuple>
namespace fs = std::filesystem;
#include <array>
#include <atomic>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
using namespace std;
using std::unordered_map;

#include "bamUtils.h"
#include "h5Writer.h"
#include "multiMap.h"
#include "samReader.h"
#include "samWriter.h"
#include "tagReadsWithGeneExon.h"
#include "umiCorrection.h"
#include "utils.h"

#include <libx/ThreadPool.hpp>
#include <libx/Timer.hpp>
#include <spdlog/spdlog.h>

enum GeneStatus
{
    Empty = 0,
    Process,
    Finish
};

template < typename T > class GenePositions
{
public:
    GenePositions() : _status(GeneStatus::Empty) {}
    ~GenePositions() {}

public:
    // Add one element
    void add(T& v)
    {
        lock_guard< mutex > lg(_mutex);
        // if (_status != GeneStatus::Empty)
        // {
        //     cerr << "add element: " << v << " but GenePositions is already finished"
        //          << endl;
        //     return;
        // }
        _data.emplace_back(v);
    }
    // Get raw data for read only
    vector< T >& get()
    {
        return _data;
    }
    // Mark it is finished
    void setStatus(GeneStatus status)
    {
        lock_guard< mutex > lg(_mutex);
        _status = status;
    }
    void updatePositions(unsigned int lessPoint, unsigned int addition)
    {
        lock_guard< mutex > lg(_mutex);
        sort(_data.begin(), _data.end());
        transform(_data.begin(), _data.end(), _data.begin(),
                  [&](T x) { return x % addition; });
        for (auto& v : _data)
        {
            if (v < lessPoint)
            {
                v += addition;
            }
        }
    }
    GeneStatus status()
    {
        lock_guard< mutex > lg(_mutex);
        return _status;
    }

private:
    vector< T > _data;
    mutex       _mutex;
    GeneStatus  _status;
};

class Bam2Gem
{
public:
    Bam2Gem(Arguments& arguments);
    ~Bam2Gem();

    bool prepare();
    int  doWork();

private:
    void checkBamFormat();

    void         I(MultiMap& mm);
    void         O(string& filename);
    void         P(TagReadsWithGeneExon* anno);
    void         dumpMetrics(string& filename, TagReadsWithGeneExon* anno);
    bool         checkPointers();
    void         plusI();
    void         plusO();
    unsigned int umicorrection(string gene);
    pair< unsigned long, unsigned long >
         prepareProcess(TagReadsWithGeneExon* anno, unsigned int start, unsigned int step);
    void manageprocess();
    // Split reads in single gene when blocked for umi correction
    // clang-format off
    typedef unordered_map< unsigned long, unordered_map< unsigned long, pair<int,bool> > > mis_map;
    typedef unordered_map< unsigned long, unordered_map< unsigned long, unsigned long > > corr_map;
    // clang-format on
    struct CorrectionData
    {
        string  genename;
        mis_map mismatch;
        exp_map expression;
    };
    void dumpSatInterData(mis_map& umiStat);
    void assignThreadNum(unsigned int totalThreadNum);

private:
    Arguments arguments;

    atomic< unsigned int > maxReads;
    vector< string >       bamFiles;

    enum ReadStatus
    {
        Empty = 0,
        Prepare,
        Ready,
        Finish
    };
    struct Record
    {
        Record()
        {
            bamRecord = createBamRecord();
        }
        Record(const Record& rhs)
        {
            if (bamRecord == nullptr)
                bamRecord = createBamRecord();
            [[maybe_unused]] const auto _ = bam_copy1(bamRecord, rhs.bamRecord);
            status                        = rhs.status.load(std::memory_order_relaxed);
        };
        Record& operator=(const Record& rhs)
        {
            [[maybe_unused]] const auto _ = bam_copy1(bamRecord, rhs.bamRecord);
            status                        = rhs.status.load(std::memory_order_relaxed);
            return *this;
        }
        ~Record()
        {
            if (bamRecord != nullptr)
            {
                destroyBamRecord(bamRecord);
                bamRecord = nullptr;
            }
        }
        BamRecord              bamRecord = nullptr;
        atomic< unsigned int > status    = 0;
    };
    vector< Record > readsQueue;

    unique_ptr< SamReader >             samReader;
    atomic< unsigned int >              pointerI, pointerO;
    array< atomic< unsigned int >, 64 > pointerP;
    atomic< bool >                      finishI, finishP, startO;

    vector< shared_ptr< GenePositions< unsigned long > > > genePos;
    map< string, vector< pair< string, unsigned int > > >  geneEnd;
    unordered_map< string, unsigned >                      geneID;
    unsigned int                                           totalGenes;
    vector< pair< string, uint32 > >                       contigs;
    vector< string >                                       geneName;

    mutex                                            prepareMutex, expMutex, blockedMutex;
    vector< queue< string > >                        preparedGenes;
    condition_variable                               cv, blockedCV;
    unordered_map< unsigned long, pair< int, int > > expMatrix;

    // Stat metrics
    atomic< unsigned long > totalReads, filteredReads, annotatedReads, uniqReads;

    // Umi correction
    UmiCorrection< unsigned long >* umiCorrection = nullptr;

    // Status of blocked
    atomic< bool > blocked;

    // Numbers of tasks running umi correction
    atomic< int > taskNum;

    FILE* satInterDataFh = nullptr;

    struct
    {
        unsigned int readThreadNum;
        unsigned int writeThreadNum;
        unsigned int annoThreadNum;
        unsigned int umiThreadNum;
    } threadNums;

    vector< CorrectionData > correctionDatas;

    atomic< uint32_t > errorNum;

    H5Writer* h5Writer = nullptr;
};
