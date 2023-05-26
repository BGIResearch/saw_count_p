/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <fstream>

#include <spdlog/spdlog.h>

#include "samReader.h"
#include "utils.h"

class SamWriter
{
public:
    SamWriter(string filename_) : filename(filename_), out(nullptr), header_(nullptr) {}
    ~SamWriter()
    {
        if (out != nullptr)
        {
            close();
        }
    }

    int close()
    {
        if (sam_idx_save(out) < 0) 
        {
            spdlog::error("Error saving index");
        }
        // if (header_ != nullptr)
        // {
        // 	bam_hdr_destroy(header_);
        // 	header_ = nullptr;
        // }
        if (out != nullptr)
        {
            hts_close(out);
            out = nullptr;
        }
        return 0;
    }

    int init(BamHeader header)
    {
        out = hts_open(filename.c_str(), "wb");
        if (!out)
        {
            string errorMsg = "Could not open " + filename;
            spdlog::error(errorMsg);
            throw std::runtime_error(ERR_CODE + "002 " + errorMsg);
        }
        header_ = header;
        int ret = sam_hdr_write(out, header_);
        if (ret != 0)
        {
            string errorMsg = "Failed to write bam header";
            spdlog::error(errorMsg);
            throw std::runtime_error(ERR_CODE + "002 " + errorMsg);
        }

        idxname = filename + ".csi";
        if (sam_idx_init(out, header_, 14, idxname.c_str()) < 0)
        {
            string errorMsg = "Failed to initialise index: " + idxname;
            spdlog::error(errorMsg);
            throw std::runtime_error(ERR_CODE + "003 " + errorMsg);
        }

        if (hts_set_opt(out, HTS_OPT_BLOCK_SIZE, _DEFAULT_HTS_BLOCK_SIZE) != 0)
        {
            string errorMsg = "Failed to set HTS_OPT_BLOCK_SIZE as "
                              + to_string(_DEFAULT_HTS_BLOCK_SIZE);
            spdlog::error(errorMsg);
            throw std::runtime_error(ERR_CODE + "003 " + errorMsg);
        }

        spdlog::debug("SamWriter init");
        return 0;
    }

    int write(BamRecord& record)
    {
        int ret = sam_write1(out, header_, record);
        if (ret < 0)
        {
            string errorMsg = "Failed to write bam record";
            spdlog::error(errorMsg);
            throw std::runtime_error(ERR_CODE + "002 " + errorMsg);
        }

        return 0;
    }

    void setThreadPool(htsThreadPool* p)
    {
        hts_set_opt(out, HTS_OPT_THREAD_POOL, p);
    }

    int setThreads(int n)
    {
        return hts_set_threads(out, n);
    }

private:
    string filename;
    string idxname;

    // A pointer to the htslib file used to access the SAM/BAM data.
    htsFile* out;

    // A htslib header data structure obtained by parsing the header of this
    // BAM.
    bam_hdr_t* header_;
};
