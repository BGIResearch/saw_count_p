/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <set>
#include <string>
#include <unordered_map>

#include "gtfRecord.h"

enum GTF_COLUMNS
{
    CHROMOSOME,
    SOURCE,
    FEATURE,
    START,
    END,
    SCORE,
    STRAND,
    FRAME,
    ATTRIBUTE
};

class AnnoReader
{
public:
    explicit AnnoReader(std::string annoFile);
    virtual ~AnnoReader();

    // Read GTF file line by line, parse data what wen want and gather by
    // genename.
    bool
    loadAnnoFile(std::unordered_map< std::string, std::vector< GTFRecord > >& gtfMap);

protected:
    struct Record
    {
        std::string geneID;
        std::string geneName;
        std::string transID;
        std::string transName;
    };
    virtual bool parseAttributes(std::string& feature, std::string& attrs,
                                 Record& record) = 0;
    std::unordered_map< std::string, std::string > splitAttrubutes(std::string& s);

protected:
    char kvDelim;  // Delimiter of key-value pairs in attribute column

private:
    std::string      annoFile;
    char             colDelim;   // Delimiter between columns
    char             attrDelim;  // Delimiter of attribute column
    static const int GTF_TAB_NUM = 9;
};

class GTFReader : public AnnoReader
{
public:
    explicit GTFReader(std::string annoFile);

protected:
    virtual bool parseAttributes(std::string& feature, std::string& attrs,
                                 Record& record);
};

class GFFReader : public AnnoReader
{
public:
    explicit GFFReader(std::string annoFile);

protected:
    virtual bool parseAttributes(std::string& feature, std::string& attrs,
                                 Record& record);

private:
    std::unordered_map< std::string, Record > recordsMap;
};

AnnoReader* createAnnoReader(std::string annoFile);