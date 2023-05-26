/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#ifdef _WIN32
#include <direct.h>
#define popen _popen
#define pclose _pclose
#else
#include <sys/sysinfo.h>
#include <unistd.h>
#endif

#include "utils.h"
#include <iostream>
#include <sstream>
#include <string.h>

#include <chrono>
#include <filesystem>
#include <string_view>
namespace fs = std::filesystem;
using namespace std;

#include <libx/File.hpp>
#include <libx/String.hpp>

unsigned int BASE2INT(unsigned char c)
{
    unsigned int res = 0;
    switch (c)
    {
    case 'A':
        res = 0;
        break;
    case 'C':
        res = 1;
        break;
    case 'G':
        res = 2;
        break;
    case 'T':
        res = 3;
        break;
    default:
        throw std::runtime_error(ERR_CODE + "003 "
                                 + "BASE2INT only support 'ACGT' as input char");
    }
    return res;
}
unsigned long seq2ul(const std::string& seq)
{
    if (seq.size() > 32)
    {
        throw std::runtime_error(ERR_CODE + "001 "
                                 + "umi maximum support 32bp length, but got "
                                 + to_string(seq.size()) + "bp");
    }
    unsigned long res = 0;
    for (auto& s : seq)
    {
        res = BASE2INT(s) + res * BASES_NUM;
    }
    return res;
}
unsigned long seq2ul(const char* seq)
{
    string tmp(seq);
    return seq2ul(tmp);
}
std::string ul2seq(unsigned long num, int len)
{
    string res(len, 'A');
    while (--len >= 0)
    {
        res.at(len) = INT2BASE[num % BASES_NUM];
        num /= BASES_NUM;
    }
    return res;
}

std::string umi2hex(unsigned long u, int len)
{
    std::stringstream ss;
    ss.width(int((len + 1) / 2));
    ss.fill('0');
    ss << std::hex << std::uppercase << u;
    return ss.str();
}
unsigned long hex2umi(std::string& h)
{
    if (h.empty())
        return 0;
    unsigned long result = stoul(h, nullptr, 16);
    return result;
}

void parseRead(BamRecord bamRecord, unsigned int& x, unsigned int& y, unsigned long& u)
{
    std::string umi = "";
    getTag(bamRecord, UR_TAG, umi);
    u = hex2umi(umi);

    int tmpx = 0, tmpy = 0;
    getTagInt(bamRecord, CX_TAG, tmpx);
    x = tmpx;
    getTagInt(bamRecord, CY_TAG, tmpy);
    y = tmpy;
}

unsigned long encodeKey(unsigned coorX, unsigned coorY, unsigned gene)
{
    unsigned long res = coorX;
    res               = (res << 22) | coorY;
    res               = (res << 20) | gene;
    return res;
}

void decodeKey(unsigned long l, unsigned* coorX, unsigned* coorY, unsigned* gene)
{
    *gene = l & 0xFFFFF;
    l >>= 20;
    *coorY = l & 0x3FFFFF;
    l >>= 22;
    *coorX = l & 0x3FFFFF;
}

int getSystemMemoryGB()
{
    // set default memory as 64GB
    int    res = 64;
    string filename("/proc/meminfo");
    auto   parseTotalMemory = [&](string& line) {
        vector< string > vec;
        libx::split(line, ' ', vec);
        try
        {
            if (vec.size() == 3)
                res = libx::to< long >(vec[1]) / 1024 / 1024;
        }
        catch (...)
        {
            res = 64;
        }
        return false;
    };
    libx::readFile(filename, parseTotalMemory);
    return res;
}

string Arguments::str()
{
    stringstream ss;
    ss << "inputBamFiles: ";
    for (auto& f : inputBamFiles)
        ss << f << " ";
    ss << std::boolalpha << "geneAnnotationFile: " << geneAnnotationFile
       << " outputBamFile: " << outputBamFile << " outputExpFile: " << outputExpFile
       << " outputMetricsFile: " << outputMetricsFile
       << " outputSatFile: " << outputSatFile << " saveLowQualReads: " << saveLowQualReads
       << " saveDupReads: " << saveDupReads << " umiPara.enable: " << umiPara.enable
       << " umiPara.minUmiNum: " << umiPara.minUmiNum
       << " umiPara.mismatch: " << umiPara.mismatch << " umiPara.len: " << umiPara.len
       << " coreNum: " << coreNum << " mapQualThre: " << mapQualThre
       << " memoryGB: " << memoryGB << " sn: " << sn << " multiMap: " << multiMap;
    return ss.str();
}

unsigned int parseResolution(string& sn)
{
    std::unordered_map< string, unsigned int > pitch(
        { { "CL1", 900 }, { "N1", 900 },   { "V3", 715 },  { "K2", 715 },
          { "S2", 715 },  { "S1", 900 },   { "F3", 715 },  { "F1", 800 },
          { "V1", 800 },  { "DP84", 715 }, { "DP8", 850 }, { "FP2", 500 },
          { "SS2", 500 }, { "FP1", 600 },  { "E1", 700 },  { "DP40", 700 },
          { "G1", 700 },  { "A", 500 },    { "B", 500 },   { "C", 500 },
          { "D", 500 },   { "U", 715 },    { "V", 715 },   { "W", 715 },
          { "X", 715 },   { "Y", 500 } });

    unsigned int result = 0;
    // check short SN
    if (sn.size() == 6)
    {
        string chip_prefix = sn.substr(0, 1);
        if (pitch.count(chip_prefix) != 0)
        {
            result = pitch[chip_prefix];
            return result;
        }
    }

    // check long SN
    string chip_prefix = sn.substr(0, 4);
    // cout<<"search filename prefix: "<<chip_prefix<<endl;
    while (!chip_prefix.empty())
    {
        if (pitch.count(chip_prefix) != 0)
        {
            result = pitch[chip_prefix];
            break;
        }
        chip_prefix.pop_back();
    }
    // cout<<"sn: "<<result<<endl;

    return result;
}