/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include <exception>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include <htslib/thread_pool.h>

#include "annotationException.h"
#include "bam2gem.h"
#include "bamUtils.h"
#include "geneBuilder.h"
#include "gtfReader.h"
#include "multiMap.h"
#include "samReader.h"
#include "samWriter.h"
#include "utils.h"

#include <libx/String.hpp>
#include <libx/System.hpp>
#include <libx/ThreadPool.hpp>
#include <libx/Timer.hpp>

#include <zlib.h>
#define cmpFile gzFile
#define cmpOpen(x) gzopen(x, "wb")
#define cmpClose(x) gzclose(x)
#define cmpFunc(x, y) gzputs(x, y)

Bam2Gem::Bam2Gem(Arguments& arguments) : arguments(arguments), errorNum(0)
{
    bamFiles = arguments.inputBamFiles;
}

Bam2Gem::~Bam2Gem()
{
    if (arguments.umiPara.enable && (umiCorrection != nullptr))
    {
        delete umiCorrection;
        umiCorrection = nullptr;
    }
    if (!arguments.outputSatFile.empty() && (satInterDataFh != nullptr))
    {
        int ret = fclose(satInterDataFh);
        if (ret != 0)
        {
            spdlog::error("error occured close file: {}, return code: {}",
                          arguments.outputSatFile, ret);
            cerr << ERR_CODE << "002 "
                 << "error occured close file: " << arguments.outputSatFile << endl;
        }
        satInterDataFh = nullptr;
    }

    if (samReader != nullptr)
        samReader->Close();

    if (h5Writer != nullptr)
    {
        delete h5Writer;
        h5Writer = nullptr;
    }
}

bool Bam2Gem::prepare()
{
    assignThreadNum(arguments.coreNum);

    finishI  = false;
    finishP  = false;
    pointerI = 0;
    pointerO = 0;

    assert(threadNums.annoThreadNum <= 64);
    for (unsigned int i = 0; i < threadNums.annoThreadNum; ++i)
    {
        pointerP[i] = i;
    }

    samReader          = SamReader::FromFile(bamFiles[0]);
    contigs            = samReader->getContigs();
    // There maybe more than one bam file, check they have same header
    if (bamFiles.size() > 1)
    {
        for (size_t i = 1; i < bamFiles.size(); ++i)
        {
            std::unique_ptr< SamReader > reader      = SamReader::FromFile(bamFiles[i]);
            auto                         tmp_contigs = reader->getContigs();
            if (tmp_contigs != contigs)
            {
                spdlog::error("Different header of bam files: {} {}", bamFiles[0],
                              bamFiles[i]);
                throw std::runtime_error(ERR_CODE + "001 "
                                         + "Different header of input bam files!");
            }
        }
    }
    spdlog::debug("Bam contigs num:{}", contigs.size());
    
    const int bytesPerRecord = 400;
    uint64_t systemMemory = arguments.memoryGB * 9.0 / 10;
    // If system memory > 1.6T, the maxReads with int type will overflow
    uint64_t estimatedReads = systemMemory * 1e9 / bytesPerRecord;
    if (estimatedReads > UINT_MAX)
        maxReads = UINT_MAX;
    else
        maxReads = estimatedReads;
    if (systemMemory <= 0 || maxReads <= 0)
    {
        std::string error = "Memory is too low: " + std::to_string(arguments.memoryGB)+"GB";
        spdlog::error(error);
        throw std::runtime_error(ERR_CODE + "001 " + error);
    }
    spdlog::debug("reads buffer size: {} buffer memory(GB): {}", maxReads, systemMemory);
    readsQueue.resize(maxReads);

    // Check if there are tags we need in bam records
    checkBamFormat();

    totalReads = filteredReads = annotatedReads = uniqReads = 0;

    if (arguments.umiPara.enable)
    {
        umiCorrection = new UmiCorrection< unsigned long >(arguments.umiPara.minUmiNum,
                                                           arguments.umiPara.mismatch,
                                                           arguments.umiPara.len);
    }

    blocked = false;

    taskNum = 0;

    if (!arguments.outputSatFile.empty())
    {
        satInterDataFh = fopen(arguments.outputSatFile.c_str(), "w");
        if (satInterDataFh == nullptr)
        {
            std::string error = "Error opening file: " + arguments.outputSatFile;
            spdlog::error(error);
            throw std::runtime_error(ERR_CODE + "001 " + error);
        }
        string header("y x geneIndex MIDIndex readCount\n");
        size_t bytes      = header.size();
        size_t writeBytes = fwrite(header.c_str(), 1, bytes, satInterDataFh);
        if (writeBytes != bytes)
        {
            spdlog::error("error occured writing sat header, need write: {} "
                            "bytes but write: {}, ferror return: {}",
                            bytes, writeBytes, ferror(satInterDataFh));
            throw std::runtime_error(ERR_CODE + "001 " + "error occured writing header of saturation file!");
        }
    }

    // Open h5 file handle
    unsigned int resolution = parseResolution(arguments.sn);
    h5Writer                = new H5Writer(arguments.outputExpFile, resolution);
    if (h5Writer == nullptr)
    {
        std::string error = "Error opening file: " + arguments.outputExpFile;
        spdlog::error(error);
        throw std::runtime_error(ERR_CODE + "001 " + error);
    }

    return true;
}

int Bam2Gem::doWork()
{
    // Load annotations.
    TagReadsWithGeneExon tagReadsWithGeneExon(arguments.geneAnnotationFile);
    if (arguments.multiMap)
        tagReadsWithGeneExon.setMultiMap();
    if (tagReadsWithGeneExon.makeOverlapDetector() != 0)
    {
        throw std::runtime_error(
            ERR_CODE + "003 "
            + "Failed makeOverlapDetector! Please check annotation file format.");
    }

    MultiMap mm(bamFiles, arguments.coreNum);
    if (arguments.multiMap && !mm.createMap(&tagReadsWithGeneExon))
    {
        throw std::runtime_error(ERR_CODE + "003 "
                                 + "Failed create map for multi-mapping reads! Please "
                                   "check bam file format or system memory");
    }

    spdlog::info("Using threads num:{}", arguments.coreNum);

    libx::Timer totalTimer("s");

    geneEnd    = tagReadsWithGeneExon.getGeneEnd();
    geneID     = tagReadsWithGeneExon.getGeneID();
    geneName   = tagReadsWithGeneExon.getGeneName();
    totalGenes = geneID.size() + 1;

    thread threadI(&Bam2Gem::I, this, std::ref(mm));
    this_thread::sleep_for(std::chrono::milliseconds(100));
    thread threadO(&Bam2Gem::O, this, std::ref(arguments.outputBamFile));
    auto   anno = &tagReadsWithGeneExon;
    P(anno);

    threadI.join();
    threadO.join();

    dumpMetrics(arguments.outputMetricsFile, anno);

    h5Writer->dump();

    spdlog::debug("memory(GB) after process: {:.3f}", libx::getSelfMemory());

    // Check is there errors?
    if (errorNum > 0)
        throw std::runtime_error(ERR_CODE + "002 " + "error occured "
                                 + to_string(errorNum) + " times");

    return 0;
}

// Check if umi exists? If true, then set length of barcode and umi
void Bam2Gem::checkBamFormat()
{
    BamRecord bamRecord = createBamRecord();
    unique_ptr< SamReader > reader = SamReader::FromFile(bamFiles[0]);
    reader->QueryAll(bamRecord);

    std::string umi = "";
    if (!getTag(bamRecord, UR_TAG, umi))
    {
        bam_destroy1(bamRecord);
        throw runtime_error(ERR_CODE + "001 "
                            + "No UR:Z tags found in bam file! Please "
                              "check the bam format.");
    }

    int cx = 0, cy = 0;
    if (!getTagInt(bamRecord, CX_TAG, cx) || !getTagInt(bamRecord, CY_TAG, cy))
    {
        bam_destroy1(bamRecord);
        throw runtime_error(ERR_CODE + "001 "
                            + "No Cx:i Cy:i tags found in bam file! Please "
                              "check the bam format.");
    }
   
    bam_destroy1(bamRecord);
    return;
}

void Bam2Gem::I(MultiMap& mm)
{
    libx::Timer timer("s"), contigTimer;

    int          hi_index, nh;
    string       qname;
    uint64       lineNum = 0;

    if (bamFiles.size() == 1)
    {
        if (samReader->setThreads(threadNums.readThreadNum) != 0)
        {
            spdlog::warn("fail init hts threads for sam reader");
        }

        auto record = createBamRecord();
        int last_contig = -1;
       
        while (true)
        {
            while (readsQueue[pointerI].status != ReadStatus::Empty)
            {
                this_thread::yield();
            }

            auto& read = readsQueue[pointerI];
            if (!samReader->QueryAll(read.bamRecord))
                break;

            // For counting total reads
            if (getTagInt(read.bamRecord, "HI", hi_index) && hi_index == 1)
                totalReads++;

            if (arguments.multiMap)
            {
                if (getTagInt(read.bamRecord, "NH", nh) && nh > 1)
                {
                    // correct read to unique from multimapping reads
                    getQName(read.bamRecord, qname);
                    int hi = mm.search(0, qname);
                    if (hi == hi_index)
                    {
                        read.bamRecord->core.qual = 255;
                        read.bamRecord->core.flag &= ~BAM_FSECONDARY;
                    }
                    else
                    {
                        read.bamRecord->core.qual = 0;
                        read.bamRecord->core.flag |= BAM_FSECONDARY;
                    }
                }
            }
            else
            {
                if (hi_index > 1)
                    continue;
            }

            read.status = ReadStatus::Prepare;
            plusI();
            lineNum++;

            if (last_contig != read.bamRecord->core.tid)
            {
                if (last_contig != -1)
                {
                    spdlog::debug("read contig: {} total reads: {} line number: {} time(s): {}",
                        contigs[last_contig].first, totalReads, lineNum, contigTimer.toc()); 
                }
                last_contig = read.bamRecord->core.tid;
            }
        }
        if (last_contig != -1)
        {
            spdlog::debug("read contig: {} total reads: {} line number: {} time(s): {}",
                contigs[last_contig].first, totalReads, lineNum, contigTimer.toc()); 
        }

        destroyBamRecord(record);
    }
    else
    {
        // Init multi samReaders and barRecords
        const int block_size = _DEFAULT_HTS_BLOCK_SIZE / bamFiles.size();
        std::vector< std::unique_ptr< SamReader > > readers;
        for (auto& filename : bamFiles)
        {
            auto reader = SamReader::FromFile(filename, block_size);
            readers.emplace_back(std::move(reader));
        }
        vector< BamRecord > records(readers.size(), nullptr);

        // TODO: threads will increase N if there are N input bam files
        // keep balance between threads resources and speed
        htsThreadPool p = { NULL, 0 };
        p.pool          = hts_tpool_init(threadNums.readThreadNum);
        for (auto& reader : readers)
        {
            if (reader->setThreadPool(&p) != 0)
            {
                spdlog::warn("fail init hts thread pool for sam reader");
            }
        }

        for (unsigned int idx = 0; idx < readers.size(); ++idx)
        {
            records[idx] = createBamRecord();
            if (!readers[idx]->QueryAll(records[idx]))
            {
                destroyBamRecord(records[idx]);
                records[idx] = nullptr;
            }
        }

        int contig_size = contigs.size();
        for (int contig = 0; contig < contig_size; ++contig)
        {
            while (true)
            {
                // 1. find maximum refStart in records
                int minStart = INT_MAX;
                int minIdx   = -1;
                for (unsigned int idx = 0; idx < records.size(); ++idx)
                {
                    auto& record = records[idx];
                    if (record == nullptr || record->core.tid != contig)
                        continue;
                    int start = getRefStart(record);
                    if (start < minStart)
                    {
                        minStart = start;
                        minIdx   = idx;
                    }
                }
                // no reads for consume
                if (minIdx == -1)
                    break;

                // 2. copy to the data queue
                auto& record = records[minIdx];
                while (readsQueue[pointerI].status != ReadStatus::Empty)
                {
                    this_thread::yield();
                }
                // For counting total reads
                if (getTagInt(record, "HI", hi_index) && hi_index == 1)
                    totalReads++;
                if (arguments.multiMap)
                {
                    auto&                       read = readsQueue[pointerI];
                    [[maybe_unused]] const auto _    = bam_copy1(read.bamRecord, record);
                    getTagInt(read.bamRecord, "NH", nh);
                    if (nh > 1)
                    {
                        // correct read to unique from multimapping reads
                        getQName(read.bamRecord, qname);
                        int hi = mm.search(minIdx, qname);
                        if (hi == hi_index)
                        {
                            read.bamRecord->core.qual = 255;
                            read.bamRecord->core.flag &= ~BAM_FSECONDARY;
                        }
                        else
                        {
                            read.bamRecord->core.qual = 0;
                            read.bamRecord->core.flag |= BAM_FSECONDARY;
                        }
                    }
                    read.status = ReadStatus::Prepare;
                    plusI();
                }
                else
                {
                    if (hi_index == 1)
                    {
                        auto&                    read = readsQueue[pointerI];
                        [[maybe_unused]] const auto _ = bam_copy1(read.bamRecord, record);
                        read.status                   = ReadStatus::Prepare;
                        plusI();
                    }
                }

                lineNum++;

                // 3. get next record
                if (!readers[minIdx]->QueryAll(record))
                {
                    destroyBamRecord(record);
                    record = nullptr;
                }
            }

            spdlog::debug("read contig: {} total reads: {} line number: {} time(s): {}",
                          contigs[contig].first, totalReads, lineNum, contigTimer.toc());
        }

        // Clear resources
        for (auto& reader : readers)
            reader->Close();
        if (p.pool)
        {
            hts_tpool_destroy(p.pool);
        }
    }

    finishI = true;
    spdlog::debug("read total reads: {} line number: {} time(s): {}", totalReads, lineNum,
                  timer.toc());
}

void Bam2Gem::O(string& filename)
{
    libx::Timer timer("s");
    unsigned    totalReads = 0;

    SamWriter samWriter(filename);
    // Update header of out bam file, adding map of gene annotation
    auto   header = samReader->getHeader();
    string co("@CO\tXF_DICT\tEXONIC:0\tINTRONIC:1\tINTERGENIC:2\n");
    char*  newt  = ( char* )realloc(header->text, header->l_text + co.size());
    header->text = newt;
    memcpy(header->text + header->l_text, co.c_str(), co.size());
    header->l_text += co.size();
    samWriter.init(header);

    if (samWriter.setThreads(threadNums.writeThreadNum) != 0)
    {
        spdlog::warn("fail init hts threads for sam writer");
    }

    auto discardLowQualReads = !arguments.saveLowQualReads;
    auto discardDupReads     = !arguments.saveDupReads;
    while (true)
    {
        // spdlog::debug("o {} {}", pointerO, readsQueue[pointerO].status);
        if (readsQueue[pointerO].status == ReadStatus::Finish)
        {
            if (discardLowQualReads && getQcFail(readsQueue[pointerO].bamRecord))
            {
            }
            else if (discardDupReads && getDuplication(readsQueue[pointerO].bamRecord))
            {
            }
            else
            {
                samWriter.write(readsQueue[pointerO].bamRecord);
            }
            readsQueue[pointerO].status = ReadStatus::Empty;
            plusO();
            totalReads++;
        }
        else
        {
            if (finishP)
                break;

            // check is blocked
            if (checkPointers() && taskNum == 0 && readsQueue[pointerO].status != 0)
            {
                // Get the chr_gene in blocked
                auto&  bamRecord = readsQueue[pointerO].bamRecord;
                string chr, ge, genename;
                chr = contigs[bamRecord->core.tid].first;
                getTag(bamRecord, GE_TAG, ge);
                // skip pseudo wakeup
                if (ge.empty())
                    continue;

                spdlog::debug("find blocked. pointer IO {} {}", pointerI, pointerO);

                blocked = true;

                genename = chr + "_" + ge;
                spdlog::debug("unblock {} read status {}", genename,
                              readsQueue[pointerO].status);

                // Do umi correction with blocked gene
                uniqReads += umicorrection(genename);

                // send signal for restart prepare thread and output thread
                // unique_lock< mutex > lk(blockedMutex);
                blocked = false;
                // blockedCV.notify_all();
            }

            this_thread::yield();
        }
    }
    samWriter.close();

    spdlog::debug("write total reads: {} time(s): {}", totalReads, timer.toc());
}
void Bam2Gem::P(TagReadsWithGeneExon* anno)
{
    libx::Timer  timer("s");
    unsigned int start = 0, step = threadNums.annoThreadNum;

    preparedGenes.resize(threadNums.annoThreadNum);
    for (unsigned int i = 0; i < totalGenes; ++i)
    {
        shared_ptr< GenePositions< unsigned long > > tmp(
            new GenePositions< unsigned long >());
        genePos.push_back(tmp);
    }

    libx::ThreadPool                                         tp(threadNums.annoThreadNum);
    vector< future< pair< unsigned long, unsigned long > > > results;
    for (; start < step; ++start)
    {
        results.emplace_back(tp.commit(
            std::bind(&Bam2Gem::prepareProcess, this, std::ref(anno), start, step)));
    }

    manageprocess();

    for (auto& result : results)
    {
        auto p = result.get();
        filteredReads += p.first;
        annotatedReads += p.second;
    }

    finishP = true;

    spdlog::debug("filterted reads: {} annotated reads: {} time(s): {}", filteredReads,
                  annotatedReads, timer.toc());
}

void Bam2Gem::dumpMetrics(string& filename, TagReadsWithGeneExon* anno)
{
    float filter_rate =
        totalReads != 0 ? float(totalReads - filteredReads) * 100 / totalReads : 0;
    float fail_annotate_rate =
        filteredReads != 0 ? float(filteredReads - annotatedReads) * 100 / filteredReads
                           : 0;
    float dup_rate = annotatedReads != 0
                         ? float(annotatedReads - uniqReads) * 100 / annotatedReads
                         : 0;
    spdlog::info("Total reads:{} Pass filter reads:{} Annotated reads:{} "
                 "Unique reads:{}",
                 totalReads, filteredReads, annotatedReads, uniqReads);
    spdlog::info("Failed filter rate:{:.2f}%", filter_rate);
    spdlog::info("Failed annotate rate:{:.2f}%", fail_annotate_rate);
    spdlog::info("Duplication rate:{:.2f}%", dup_rate);

    // Dump metrics information to disk.
    std::string   anno_metrics = anno->dumpMetrics();
    std::ofstream ofs(filename, std::ofstream::out);
    ofs.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    if (ofs.is_open())
    {
        try
        {
            ofs.precision(2);
            ofs.setf(std::ios::fixed);
            // clang-format off
            ofs << "## FILTER & DEDUPLICATION METRICS\n"
                "TOTAL_READS\tPASS_FILTER\tANNOTATED_READS\tUNIQUE_READS\t"
                "FAIL_FILTER_RATE\tFAIL_ANNOTATE_RATE\tDUPLICATION_RATE\n";
            // clang-format on
            ofs << totalReads << "\t" << filteredReads << "\t" << annotatedReads << "\t"
                << uniqReads << "\t" << filter_rate << "\t" << fail_annotate_rate << "\t"
                << dup_rate << std::endl;
            ofs << anno_metrics;
            ofs.close();
            spdlog::info("Success dump metrics file: {}", filename);
        }
        catch (std::ofstream::failure& e)
        {
            string errorMsg = "Error writing file: " + filename;
            spdlog::error(errorMsg);
            throw std::runtime_error(ERR_CODE + "002 " + errorMsg);
        }
    }
    else
    {
        string errorMsg = "Error opening file: " + filename;
        spdlog::error(errorMsg);
        throw std::runtime_error(ERR_CODE + "002 " + errorMsg);
    }
}

// If there is blocked, the process thread pointers must
// in range [pointerI,pointerI+1,...,pointerI+N]
bool Bam2Gem::checkPointers()
{
    if (pointerI != pointerO)
        return false;

    set< unsigned int > pointerSet;
    for (unsigned int i = 0; i < threadNums.annoThreadNum; ++i)
        pointerSet.insert((pointerO + i) % maxReads);

    for (unsigned int i = 0; i < threadNums.annoThreadNum; ++i)
    {
        unsigned int p = pointerP[i];
        // spdlog::debug("p: {} {}", i, p);
        if (pointerSet.count(p) == 0)
        {
            return false;
        }
    }
    return true;
}

void Bam2Gem::plusI()
{
    if (++pointerI == maxReads)
        pointerI = 0;
}
void Bam2Gem::plusO()
{
    if (++pointerO == maxReads)
        pointerO = 0;
}
unsigned int Bam2Gem::umicorrection(string gene)
{
    // spdlog::debug("begin umicorrection gene: {}", gene);

    unsigned int result = 0;
    // vector<unsigned int> readslist;

    unsigned int geneid = geneID[gene];
    // for (auto& partPos : genePos)
    // {
    //     if (partPos[geneid].empty()) continue;
    //     auto& partlist = partPos[geneid];
    //     readslist.insert(readslist.end(), partlist.begin(), partlist.end());
    //     partlist.clear();
    //     partlist.shrink_to_fit();
    // }
    // if (readslist.empty()) return result;

    genePos[geneid]->setStatus(GeneStatus::Process);

    auto& temp = genePos[geneid]->get();

    // if (temp.empty())
    //     return result;

    taskNum += 1;
    remove_reference< decltype(temp) >::type readslist;
    readslist.swap(temp);
    std::sort(readslist.begin(), readslist.end());
    std::transform(readslist.begin(), readslist.end(), readslist.begin(),
                   [&](unsigned long x) { return x % maxReads; });

    mis_map umi_mismatch;
    exp_map barcode_gene_exp;
    // Check is there unfinished correction data in memory
    for (auto& corrData : correctionDatas)
    {
        if (corrData.genename == gene)
        {
            spdlog::debug("gene {} retrieve correction data", gene);
            umi_mismatch.swap(corrData.mismatch);
            barcode_gene_exp.swap(corrData.expression);
            break;
        }
    }

    unsigned int  coorX, coorY;
    unsigned long umi;
    for (auto& pos : readslist)
    {
        auto& bamRecord = readsQueue[pos].bamRecord;

        parseRead(bamRecord, coorX, coorY, umi);

        int locus;
        getTagInt(bamRecord, FUNCTION_TAG, locus);
        bool exonic = (locus == 0);
        auto key    = encodeKey(coorX, coorY, geneid);
        if (umi_mismatch.count(key) == 0)
        {
            umi_mismatch[key] = {};
        }
        if (umi_mismatch[key].count(umi) == 0)
        {
            umi_mismatch[key][umi] = { 1, exonic };
        }
        else
        {
            umi_mismatch[key][umi].first++;
            umi_mismatch[key][umi].second &= exonic;
            // Save the reads that duplicate
            setDuplication(bamRecord);
        }

        result++;
    }

    corr_map umi_correct;
    if (arguments.umiPara.enable)
        umiCorrection->deDupUmi(umi_mismatch, umi_correct);
    auto umiLen = arguments.umiPara.len;

    for (auto& pos : readslist)
    {
        auto& bamRecord = readsQueue[pos].bamRecord;
        if (getDuplication(bamRecord))
        {
            result--;
            readsQueue[pos].status = ReadStatus::Finish;
            continue;
        }

        parseRead(bamRecord, coorX, coorY, umi);
        auto key = encodeKey(coorX, coorY, geneid);
        if (umi_mismatch[key][umi].first == 0)
        {
            result--;
            [[maybe_unused]] int ret     = 0;
            auto                 seq     = umi_correct[key][umi];
            std::string          correct = umi2hex(seq, umiLen);
            ret = bam_aux_append(bamRecord, UB_TAG, 'Z', correct.size() + 1,
                                 ( uint8_t* )correct.c_str());
            setDuplication(bamRecord);
        }
        else
        {
            string strand;
            getTag(bamRecord, STRAND_TAG, strand);
            bool geneNegative = strand == "-";
            if (geneNegative == getNegativeStrand(bamRecord))
            {
                // Calculate barcode gene expression
                barcode_gene_exp[key].first++;
                if (umi_mismatch[key][umi].second)
                    barcode_gene_exp[key].second++;
            }
        }
        readsQueue[pos].status = ReadStatus::Finish;
    }

    lock_guard< mutex > lck(expMutex);
    if (blocked)
    {
        // blocked, save correction data to memory
        // find position to store data
        spdlog::debug("gene {} store correction data in memory", gene);
        unsigned int pos = 0;
        for (; pos < correctionDatas.size(); ++pos)
        {
            if (correctionDatas[pos].genename == gene)
                break;
        }
        if (correctionDatas.empty() || pos == correctionDatas.size())
        {
            // add a new element
            correctionDatas.push_back(CorrectionData());
        }
        auto& corrData    = correctionDatas[pos];
        corrData.genename = gene;
        corrData.mismatch.swap(umi_mismatch);
        corrData.expression.swap(barcode_gene_exp);
    }
    else
    {
        // not blocked, save correction data to file
        if (barcode_gene_exp.size() > 0)
            h5Writer->store(geneName[geneid], barcode_gene_exp);

        if (!arguments.outputSatFile.empty() && (satInterDataFh != nullptr)
            && umi_mismatch.size() > 0)
            dumpSatInterData(umi_mismatch);

        genePos[geneid]->setStatus(GeneStatus::Finish);
    }

    taskNum -= 1;
    // spdlog::debug("end umicorrection gene: {}", gene);
    return result;
}

pair< unsigned long, unsigned long >
Bam2Gem::prepareProcess(TagReadsWithGeneExon* anno, unsigned int start, unsigned int step)
{
    spdlog::debug("start anno thread: {}", start);
    libx::Timer                            timer("s");
    unsigned                               totalReads = 0;
    vector< pair< string, unsigned int > > ends;
    string                                 currCtg;
    auto&                                  partPos     = genePos;
    auto&                                  partGenes   = preparedGenes[start];
    unsigned int                           notifyTimes = 0;
    unsigned long                          addition    = 0;
    unsigned long                          filtered = 0, annotated = 0;
    auto&                                  pos         = pointerP[start];
    auto                                   mapQualThre = arguments.mapQualThre;

    while (true)
    {
        if (readsQueue[pos].status == ReadStatus::Prepare)
        {
            auto& bamRecord = readsQueue[pos].bamRecord;

            string ctg = contigs[bamRecord->core.tid].first;
            if (currCtg.empty() || currCtg != ctg)
            {
                while (ends.size() > 0)
                {
                    string genename = ends.front().first;
                    ends.erase(ends.begin());

                    unique_lock< mutex > lk(prepareMutex);
                    partGenes.push(genename);
                    cv.notify_one();
                    notifyTimes++;
                }
                currCtg = ctg;
                // incase there are contigs in bam but not in gtf/gff
                if (geneEnd.count(currCtg) != 0)
                    ends = geneEnd[currCtg];
            }

            // Filter mapping quality.
            int score = getQual(bamRecord);
            if (score < mapQualThre)
            {
                setQcFail(bamRecord);
                readsQueue[pos].status = ReadStatus::Finish;
            }
            else
            {
                filtered++;

                unsigned int strand = 0, locus = 0;
                string       ge_value;
                // Set annotations, need the gene name for the next step
                anno->setAnnotation(bamRecord, ctg, ge_value, &strand, &locus);
                if (locus != 0)
                {
                    [[maybe_unused]] int ret =
                        updateIntTags(bamRecord, FUNCTION_TAG, LocusTransfer[locus]);
                }
                if (strand != 0)
                {
                    annotated++;

                    updateStrTags(bamRecord, GE_TAG, ge_value);
                    updateStrTags(bamRecord, STRAND_TAG, strand == 1 ? "+" : "-");
                    readsQueue[pos].status = ReadStatus::Ready;
                    unsigned int  geneid   = geneID[currCtg + "_" + ge_value];
                    unsigned long realPos  = pos + addition;
                    partPos[geneid]->add(realPos);
                    // Start a new thread for correcting umi
                    while (
                        (ends.size() > 0)
                        && (( unsigned int )getRefStart(bamRecord) > ends.front().second))
                    {
                        string genename = ends.front().first;
                        ends.erase(ends.begin());

                        unique_lock< mutex > lk(prepareMutex);
                        partGenes.push(genename);
                        cv.notify_one();
                        notifyTimes++;
                    }
                }
                else
                    readsQueue[pos].status = ReadStatus::Finish;
            }

            pos += step;
            // spdlog::debug("start anno thread: {} {}", start, pos);
            if (pos >= maxReads)
            {
                pos -= maxReads;
                addition += maxReads;
            }
            totalReads++;
        }
        else
        {
            if (finishI)
                break;
            // this_thread::yield();
            // if (multiple == 0)
            //     this_thread::sleep_for(std::chrono::milliseconds(100));
            // else
            //     this_thread::yield();
            this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    while (ends.size() > 0)
    {
        string genename = ends.front().first;
        ends.erase(ends.begin());

        unique_lock< mutex > lk(prepareMutex);
        partGenes.push(genename);
        cv.notify_one();
        notifyTimes++;
    }

    unique_lock< mutex > lk(prepareMutex);
    partGenes.push("");
    cv.notify_one();

    spdlog::debug("prepare total reads: {} time(s): {}", totalReads, timer.toc());
    return { filtered, annotated };
}
void Bam2Gem::manageprocess()
{
    libx::Timer timer("s");

    queue< string > orderedGenes;

    for (auto& [contig, _] : contigs)
    {
        if (geneEnd.count(contig) == 0)
            continue;

        for (auto& [genename, _] : geneEnd[contig])
            orderedGenes.push(genename);
    }

    libx::ThreadPool                 tp(threadNums.umiThreadNum);
    vector< future< unsigned int > > results;
    while (true)
    {
        unique_lock< mutex > lk(prepareMutex);

        cv.wait(lk, [&] {
            for (auto& partQueue : preparedGenes)
            {
                if (partQueue.empty())
                    return false;
            }
            return true;
        });

        string genename;
        while (!orderedGenes.empty())
        {
            string currGene = orderedGenes.front();
            for (auto& partQueue : preparedGenes)
            {
                if (partQueue.front() == currGene)
                {
                    partQueue.pop();
                    genename = currGene;
                }
            }
            if (!genename.empty())
                break;
            orderedGenes.pop();
        }

        // End prepare threads
        if (genename.empty())
            break;

        results.emplace_back(
            tp.commit(std::bind(&Bam2Gem::umicorrection, this, genename)));
    }

    for (auto& result : results)
    {
        uniqReads += result.get();
    }
    spdlog::debug("uniq reads: {} time(s): {}", uniqReads, timer.toc());
}

void Bam2Gem::dumpSatInterData(mis_map& umiStat)
{
    // Just store data to file
    char         buf[1024];
    unsigned int gene, coorX, coorY;
    for (const auto& [b, p] : umiStat)
    {
        decodeKey(b, &coorX, &coorY, &gene);
        for (const auto& [umi, count] : p)
        {
            if (count.first == 0)
                continue;
            sprintf(buf, "%d %d %d %lu %d\n", coorY, coorX, gene, umi, count.first);
            size_t bytes      = strlen(buf);
            size_t writeBytes = fwrite(buf, 1, bytes, satInterDataFh);
            if (writeBytes != bytes)
            {
                spdlog::error("error occured writing sat inter data, need write: {} "
                              "bytes but write: {}, ferror return: {}",
                              bytes, writeBytes, ferror(satInterDataFh));
                errorNum += 1;
                return;
            }
        }
    }
}

void Bam2Gem::assignThreadNum([[maybe_unused]] unsigned int totalThreadNum)
{
    // there are fixed thread numbers: read/write/read daemon
    unsigned int coef = (totalThreadNum - 3) / 5;
    if (coef < 1)
    {
        string errorMsg = "Too less threads for runing, at least 8 but got "
                          + to_string(totalThreadNum);
        spdlog::warn(errorMsg);
        coef = 1;
        // throw std::runtime_error(ERR_CODE+"001 "+errorMsg);
    }
    threadNums.readThreadNum  = 1 * coef;
    threadNums.writeThreadNum = 2 * coef;
    threadNums.annoThreadNum  = 1 * coef;
    threadNums.umiThreadNum   = 1 * coef;

    spdlog::debug("assign threads, read bam: {} write bam: {} annotation: {} "
                  "umi correction: {}",
                  threadNums.readThreadNum + 2, threadNums.writeThreadNum + 2,
                  threadNums.annoThreadNum, threadNums.umiThreadNum);
}
