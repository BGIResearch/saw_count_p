/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#pragma once

#include "utils.h"

#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

template < typename T > class UmiCorrection
{
public:
    typedef T umi_type;

    UmiCorrection() : UmiCorrection(5, 1, 10) {}
    UmiCorrection(int numThre, int mismatch, int umiLen)
        : numThre(numThre), mismatch(mismatch), umiLen(umiLen)
    {
    }
    ~UmiCorrection() {}

    int
    deDupUmi(std::unordered_map< unsigned long,
                                 std::unordered_map< umi_type, pair< int, bool > > >&
                                                                             umi_mismatch,
             std::unordered_map< unsigned long,
                                 std::unordered_map< umi_type, umi_type > >& umi_correct)
    {
        if (umi_mismatch.empty())
            return 0;

        vector< pair< umi_type, int > > array;
        for (auto& p : umi_mismatch)
        {
            if (p.second.size() < numThre)
                continue;

            array.clear();
            // Transform data from map to vector<pair> for sorting by value
            for (const auto& umi : p.second)
            {
                array.push_back({ umi.first, umi.second.first });
            }
            // Sort the vector by cnt
            sort(array.begin(), array.end(), compareBySecond);
            // Pairwise comparison of all umis
            for (size_t i = array.size() - 1; i > 0; --i)
            {
                for (size_t j = 0; j < i; ++j)
                {
                    // Calculate distance of two umis
                    auto umi1 = array[i].first;
                    auto umi2 = array[j].first;
                    if (similar(umi1, umi2))
                    {
                        p.second[umi2].first += p.second[umi1].first;
                        p.second[umi2].second &= p.second[umi1].second;
                        p.second[umi1].first = 0;

                        // Mark the correct umi
                        if (umi_correct.count(p.first) == 0)
                            umi_correct[p.first] = {};
                        umi_correct[p.first][umi1] = umi2;
                        for (auto& [k, v] : umi_correct[p.first])
                        {
                            if (v == umi1)
                                v = umi2;
                        }

                        // Break current umi
                        break;
                    }
                }
            }
        }

        return 0;
    }
    // Calculate mismatch of two umis
    // Return true: used for correction, false: not similar
    bool similar(umi_type v1, umi_type v2)
    {
        unsigned int distance = 0;
        for (size_t i = 0; i < umiLen; ++i)
        {
            if ((v1 & 0x3) != (v2 & 0x3))
            {
                if ((++distance) > mismatch)
                    return false;
            }
            v1 >>= 2;
            v2 >>= 2;
        }
        return true;
    }

private:
    inline static bool compareBySecond(const pair< umi_type, int >& p1,
                                       const pair< umi_type, int >& p2)
    {
        if (p1.second > p2.second)
            return true;
        else if (p1.second < p2.second)
            return false;
        else
            return p1.first > p2.first;
    }

private:
    unsigned int numThre;
    unsigned int mismatch;
    unsigned int umiLen;
};