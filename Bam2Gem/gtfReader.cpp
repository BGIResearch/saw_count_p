/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include <filesystem>
namespace fs = std::filesystem;
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#include <libx/File.hpp>
#include <libx/GZFile.hpp>
#include <libx/String.hpp>
#include <spdlog/spdlog.h>

#include "annotationException.h"
#include "gtfReader.h"

AnnoReader::AnnoReader(string annoFile) : annoFile(annoFile)
{
    colDelim  = '\t';
    attrDelim = ';';
}

AnnoReader::~AnnoReader() {}

bool AnnoReader::loadAnnoFile(unordered_map< string, vector< GTFRecord > >& gtfMap)
{
    // Travel the records by line.
    int              num = 0;
    vector< string > splitedLine;
    auto             processLine = [&](string& line) {
        if (line.empty() || line[0] == '#')
            return true;

        // Parse data and construct instance of GTFRecord.
        // Format:CHROMOSOME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME,
        // ATTRIBUTE
        splitedLine.clear();
        libx::split(line, colDelim, splitedLine);
        if (splitedLine.size() != GTF_TAB_NUM)
            return true;

        // Check if the record is valid
        Record record;
        if (!parseAttributes(splitedLine[GTF_COLUMNS::FEATURE],
                             splitedLine[GTF_COLUMNS::ATTRIBUTE], record))
            return true;

        GTFRecord gtfRecord(
            splitedLine[GTF_COLUMNS::CHROMOSOME], splitedLine[GTF_COLUMNS::FEATURE],
            stoi(splitedLine[GTF_COLUMNS::START]), stoi(splitedLine[GTF_COLUMNS::END]),
            splitedLine[GTF_COLUMNS::STRAND] == "-", record.geneID, record.geneName,
            record.transID, record.transName);
        // Gather gtfRecord by genename.
        // Skip pseudo genes with no gene name
        if (++num % 100000 == 0)
        {
            spdlog::get("gtf")->info(
                "read {:10d} valid GTF records. Last read position: {}:{}", num,
                gtfRecord.getContig(), gtfRecord.getStart());
        }
        if (!gtfRecord.getGeneName().empty())
        {
            gtfMap[gtfRecord.getGeneName()].push_back(std::move(gtfRecord));
        }

        return true;
    };

    if (libx::endswith(annoFile, ".gz"))
        libx::readGZFile(annoFile, processLine);
    else
        libx::readFile(annoFile, processLine);

    return true;
}

unordered_map< string, string > AnnoReader::splitAttrubutes(string& s)
{
    unordered_map< string, string > res;
    vector< string >                splited;
    libx::split(s, attrDelim, splited);
    for (auto& p : splited)
    {
        p = libx::trimWhitespace(p);
        string k, v;
        libx::split(p, kvDelim, k, v);
        k      = libx::trimWhitespace(k);
        v      = libx::trimWhitespace(v);
        v      = libx::trim(v, [](const char c) { return (c == '"'); });
        res[k] = v;
    }
    return res;
}

GTFReader::GTFReader(string annoFile) : AnnoReader(annoFile)
{
    kvDelim = ' ';
}

bool GTFReader::parseAttributes(std::string& feature, std::string& attrs, Record& res)
{
    auto kvpairs = splitAttrubutes(attrs);

    string id, name;
    id   = kvpairs.count("gene_id") ? kvpairs.at("gene_id") : "";
    name = kvpairs.count("gene_name") ? kvpairs.at("gene_name") : "";
    if (id.empty() || name.empty())
        return false;

    res.geneID   = id;
    res.geneName = name;
    if (feature != "gene")
    {
        id   = kvpairs.count("transcript_id") ? kvpairs.at("transcript_id") : "";
        name = kvpairs.count("transcript_name") ? kvpairs.at("transcript_name") : "";
        if (id.empty() || name.empty())
            return false;

        res.transID   = id;
        res.transName = name;
    }

    return true;
}

GFFReader::GFFReader(string annoFile) : AnnoReader(annoFile)
{
    kvDelim = '=';
}

bool GFFReader::parseAttributes(string& feature, string& attrs, Record& res)
{
    auto   kvpairs = splitAttrubutes(attrs);
    string id, name, parent;
    id     = kvpairs.count("ID") ? kvpairs.at("ID") : "";
    name   = kvpairs.count("Name") ? kvpairs.at("Name") : "";
    parent = kvpairs.count("Parent") ? kvpairs.at("Parent") : "";
    if ((id.empty() || name.empty()) && feature != "exon")
    {
        // spdlog::debug("skip non-exon line which lack ID or Name: {}", attrs);
        return false;
    }

    if (feature == "gene")
    {
        res.geneID     = id;
        res.geneName   = name;
        recordsMap[id] = res;  // Store data for searching parent
    }
    else
    {
        vector< string > tempParents;
        libx::split(parent, ',', tempParents);
        for (auto& tp : tempParents)
        {
            // There is no parent of this transcript
            if (!recordsMap.count(tp))
                continue;

            auto& record = recordsMap.at(tp);
            res.geneID   = record.geneID;
            res.geneName = record.geneName;
            // The parent is gene
            if (feature == "mRNA")
            {
                res.transID   = id;
                res.transName = name;
            }
            else
            {
                res.transID   = record.transID;
                res.transName = record.transName;
            }
        }
        if (!res.geneID.empty() && !res.geneName.empty() && !res.transID.empty()
            && !res.transName.empty())
            recordsMap[id] = res;
        else
            return false;
    }
    return true;
}

AnnoReader* createAnnoReader(string annoFile)
{
    if (libx::endswith(annoFile, ".gtf") || libx::endswith(annoFile, ".gtf.gz"))
        return new GTFReader(annoFile);
    else if (libx::endswith(annoFile, ".gff") || libx::endswith(annoFile, ".gff3")
             || libx::endswith(annoFile, ".gff.gz")
             || libx::endswith(annoFile, ".gff3.gz"))
        return new GFFReader(annoFile);
    else
        throw AnnotationException("Invalid gtf/gff file extension: " + annoFile
                                  + ", support .gtf/.gtf.gz/.gff/.gff.gz/.gff3/.gff3.gz");
}