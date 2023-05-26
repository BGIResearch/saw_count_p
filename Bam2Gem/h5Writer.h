/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#pragma once

#include "hdf5.h"
#include <spdlog/spdlog.h>

#include "utils.h"

#include <list>
#include <map>
#include <string>
#include <vector>

class H5Writer
{
public:
    H5Writer(std::string filename, unsigned int resolution);
    ~H5Writer();

    // Store geneExp->bin1->expression by gene
    bool store(std::string& genename, exp_map& em);

    // Dump to disk
    bool dump();

private:
    static constexpr size_t char_len = 32;
    struct Dnb
    {
        Dnb(unsigned int x, unsigned int y, unsigned int cnt) : x(x), y(y), cnt(cnt) {}
        unsigned int x;
        unsigned int y;
        unsigned int cnt;
    };
    struct GenePos
    {
        GenePos(const char* g, unsigned int o, unsigned c)
        {
            int i = 0;
            while (g[i] != '\0')
            {
                gene[i] = g[i];
                ++i;
            }
            offset = o;
            count  = c;
        }
        char         gene[char_len] = { 0 };
        unsigned int offset;
        unsigned int count;
    };

    std::vector< Dnb >           dnbList;
    std::vector< GenePos >       genePosList;
    std::map< std::string, int > genes;
    std::vector< unsigned int >  exonList;

    // ids of hdf5 library
    hid_t  m_file_id;
    hid_t  m_group_id;
    hid_t  m_dataspace_id;
    hid_t  m_dataset_id;
    herr_t m_status;

    // stat values
    int          maxexp, maxexon;
    unsigned int min_x, min_y, max_x, max_y;
    unsigned int resolution;
};