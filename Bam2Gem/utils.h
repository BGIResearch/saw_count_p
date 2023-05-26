/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <sys/syscall.h>
#define gettid() syscall(__NR_gettid)

#include <exception>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

#include "bamUtils.h"

typedef signed char int8;
typedef short       int16;
typedef int         int32;
typedef long long   int64;

typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;

typedef unordered_map< unsigned long, std::pair< int, int > > exp_map;

const int _DEFAULT_HTS_BLOCK_SIZE = 128 * (1024 * 1024);

static const char HI_TAG[]       = "HI";  // Hit index
static const char UR_TAG[]       = "UR";  // Umi sequence
static const char UB_TAG[]       = "UB";  // Umi sequence after correction
static const char GE_TAG[]       = "GE";  // Gene name
static const char FUNCTION_TAG[] = "XF";  // Locus function
static const char STRAND_TAG[]   = "GS";  // Strand
static const char CX_TAG[]       = "Cx";  // Coordinate X, component of Barcode
static const char CY_TAG[]       = "Cy";  // Coordinate Y, component of Barcode

static const string ERR_CODE = "SAW-A20";

// Transform between (umi)sequence and unsigned long
static const int  BASES_NUM  = 4;
static const char INT2BASE[] = "ACGT";
unsigned int      BASE2INT(unsigned char c);
unsigned long     seq2ul(const std::string& seq);
unsigned long     seq2ul(const char* seq);
std::string       ul2seq(unsigned long num, int len = 10);
// Parse coorX/coorY/umi from bam record according tags
void parseRead(BamRecord bamRecord, unsigned int& x, unsigned int& y, unsigned long& u);
// Transform between coorX/coorY/gene and unsigned long
unsigned long encodeKey(unsigned coorX, unsigned coorY, unsigned gene);
void decodeKey(unsigned long l, unsigned* coorX, unsigned* coorY, unsigned* gene);
// Transform between umi sequence and hex
std::string   umi2hex(unsigned long u, int len = 10);
unsigned long hex2umi(std::string& h);
// System functions
int getSystemMemoryGB();
// Parse resolution from input bam filename
// Used for output gef file
unsigned int parseResolution(string& filename);

// Parameters of user input
struct UmiParameter
{
    bool enable;  // control of umi correction, true for enable and false for
                  // disable

    int minUmiNum;  // minimum umi types for correction
    int mismatch;   // number of bases allowed to be wrong
    int len;        // length of umi sequence
};
struct Arguments
{
    vector< string > inputBamFiles;
    string           geneAnnotationFile;

    string outputBamFile;
    string outputExpFile;
    string outputMetricsFile;
    string outputSatFile;

    bool saveLowQualReads;
    bool saveDupReads;

    UmiParameter umiPara;

    int coreNum;
    int memoryGB;

    int mapQualThre;

    string sn;

    bool multiMap;

    // Serialization
    string str();
};
