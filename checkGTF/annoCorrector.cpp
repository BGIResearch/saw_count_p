/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include "annoCorrector.h"

#include <climits>

#include <libx/File.hpp>
#include <libx/GZFile.hpp>
#include <libx/String.hpp>

AnnoCorrector::AnnoCorrector() {}

AnnoCorrector::~AnnoCorrector() {}
bool AnnoCorrector::load(string filename)
{
    auto parser = [&](string& line) {
        rawRecords.emplace_back(line);
        return true;
    };
    if (libx::endswith(filename, ".gz"))
        libx::readGZFile(filename, parser, "");
    else
        libx::readFile(filename, parser, "");

    return true;
}

bool AnnoCorrector::dump(string filename)
{
    if (libx::endswith(filename, ".gz"))
        libx::writeGZFile(corrRecords, filename);
    else
        libx::writeFile(corrRecords, filename);
    return true;
}

void AnnoCorrector::checkStrand(string& strand)
{
    if (strand == "_")
        strand = "-";
}

string AnnoCorrector::getAttr(string& attrStr, string key)
{
    vector< string > attrs;
    libx::split(attrStr, ';', attrs);
    map< string, string > attrPairs;
    for (auto& p : attrs)
    {
        string key, value;
        p = libx::trimWhitespace(p);
        libx::split(p, kvDelim, key, value);
        attrPairs.insert({ key, value });
    }

    if (attrPairs.count(key))
        return attrPairs[key];
    else
        return "";
}

string AnnoCorrector::fillupAttr(string& attrStr)
{
    vector< string > attrs;
    libx::split(attrStr, ';', attrs);
    map< string, string > attrPairs;
    for (auto& p : attrs)
    {
        string key, value;
        p = libx::trimWhitespace(p);
        libx::split(p, kvDelim, key, value);
        attrPairs.insert({ key, value });
    }
    for (auto& [k1, k2] : keyMap)
    {
        if (attrPairs.count(k1) != 0 && attrPairs.count(k2) == 0)
        {
            attrPairs.insert({ k2, attrPairs[k1] });
            if (!libx::endswith(attrStr, ";") && !libx::endswith(attrStr, "; "))
                attrStr.push_back(';');
            attrStr += libx::join({ k2, attrPairs[k1] }, kvDelim);
        }
    }

    if (attrPairs.count(gidName) == 0)
        // throw std::runtime_error("not exists gene id: "+attrStr);
        return "";

    return attrPairs[gidName];
}

GTFCorrector::GTFCorrector() : AnnoCorrector()
{
    kvDelim = ' ';
    keyMap.insert({ GID, GNAME });
    keyMap.insert({ GNAME, GID });
    keyMap.insert({ TID, TNAME });
    keyMap.insert({ TNAME, TID });
    gidName = GID;
}

bool GTFCorrector::correction()
{
    vector< vector< string > > cache;
    string                     lastGid;
    for (auto& line : rawRecords)
    {
        if (line.at(0) == '#')
        {
            cache.push_back({ line });
            continue;
        }
        vector< string > sline;
        libx::split(line, '\t', sline);
        if (sline.size() != 9)
        {
            throw std::runtime_error("invalid column number: " + line);
            continue;
        }
        // correct _ to -
        checkStrand(sline[6]);
        // fill up id/name
        auto gid = fillupAttr(sline.back());
        // skip the line without id and name
        if (gid.empty())
            continue;
        // fill up lacking lines
        if (gid == lastGid || lastGid.empty())
        {
            cache.emplace_back(sline);
            lastGid = gid;
        }
        else
        {
            auto lines = fillupLine(cache);
            corrRecords.insert(corrRecords.end(), lines.begin(), lines.end());
            cache.clear();

            cache.emplace_back(sline);
            lastGid = gid;
        }
    }
    auto lines = fillupLine(cache);
    corrRecords.insert(corrRecords.end(), lines.begin(), lines.end());
    cache.clear();

    return true;
}

vector< string > GTFCorrector::fillupLine(vector< vector< string > >& raw)
{
    if (raw.empty())
        return {};

    vector< string > res;
    // check if lacking gene line or transcript line
    bool lackGeneLine = true, lackTransLine = true, lackExonLine = true;
    for (auto& line : raw)
    {
        if (line[0][0] == '#')
            continue;
        if (line[2] == "gene")
            lackGeneLine = false;
        else if (line[2] == "transcript")
            lackTransLine = false;
        else if (line[2] == "exon")
            lackExonLine = false;
    }
    if ((!lackGeneLine && !lackTransLine) || lackExonLine)
    {
        for (auto& line : raw)
        {
            if (line[0][0] == '#')
                res.push_back(line[0]);
            else
                res.push_back(libx::join(line, '\t'));
        }
    }
    else
    {
        int                        firstLinePos = 0;
        string                     gid;
        string                     gname;
        vector< vector< string > > cache;
        for (auto& line : raw)
        {
            if (line[0][0] != '#')
            {
                gid   = getAttr(line.back(), GID);
                gname = getAttr(line.back(), GNAME);
                break;
            }
            firstLinePos++;
        }

        int    geneStart = INT_MAX, geneEnd = INT_MIN;
        int    transStart = INT_MAX, transEnd = INT_MIN;
        string transID;
        int    lastTransPos = -1;
        for (auto& line : raw)
        {
            if (line[0][0] == '#')
            {
                cache.push_back(line);
                continue;
            }
            if (line[2] == "gene")
            {
                cache.push_back(line);
                continue;
            }

            string tid = getAttr(line.back(), TID);
            if (lackTransLine && (transID != tid || transID.empty()))
            {
                vector< string > transLine = line;
                transLine[2]               = "transcript";
                string tname               = getAttr(line.back(), TNAME);
                transLine.back()           = libx::join({ GID, gid }, kvDelim) + ";"
                                   + libx::join({ GNAME, gname }, kvDelim) + ";"
                                   + libx::join({ TID, tid }, kvDelim) + ";"
                                   + libx::join({ TNAME, tname }, kvDelim);
                cache.push_back(transLine);
                if (lastTransPos > 0)
                {
                    cache[lastTransPos][3] = std::to_string(transStart);
                    cache[lastTransPos][4] = std::to_string(transEnd);
                }
                lastTransPos = cache.size() - 1;
                transStart   = INT_MAX;
                transEnd     = INT_MIN;

                transID = tid;
            }
            int start  = stoi(line[3]);
            int end    = stoi(line[4]);
            geneStart  = min(geneStart, start);
            geneEnd    = max(geneEnd, end);
            transStart = min(transStart, start);
            transEnd   = max(transEnd, end);

            cache.push_back(line);
        }
        if (lastTransPos > 0)
        {
            cache[lastTransPos][3] = std::to_string(transStart);
            cache[lastTransPos][4] = std::to_string(transEnd);
        }
        if (lackGeneLine)
        {
            vector< string > geneLine = raw[firstLinePos];

            geneLine[2]     = "gene";
            geneLine[3]     = std::to_string(geneStart);
            geneLine[4]     = std::to_string(geneEnd);
            geneLine.back() = libx::join({ GID, gid }, kvDelim) + ";"
                              + libx::join({ GNAME, gname }, kvDelim);
            cache.insert(cache.begin() + firstLinePos, geneLine);
        }

        for (auto& line : cache)
        {
            if (line[0][0] == '#')
                res.push_back(line[0]);
            else
                res.push_back(libx::join(line, '\t'));
        }
    }

    return res;
}

GFFCorrector::GFFCorrector() : AnnoCorrector()
{
    kvDelim = '=';
    keyMap.insert({ ID, NAME });
    keyMap.insert({ NAME, ID });
    gidName = ID;
}

bool GFFCorrector::correction()
{
    string gid, tid;
    for (auto& line : rawRecords)
    {
        if (line.at(0) == '#')
        {
            corrRecords.push_back(line);
            continue;
        }
        vector< string > sline;
        libx::split(line, '\t', sline);
        if (sline.size() != 9)
        {
            throw std::runtime_error("invalid column number: " + line);
            continue;
        }
        // correct _ to -
        checkStrand(sline[6]);
        // fill up id/name
        auto id = fillupAttr(sline.back());
        // skip the line without ID and Name
        if (id.empty())
            continue;
        if (sline[2] == "gene")
            gid = id;
        else if (sline[2] == "mRNA")
            tid = id;
        // fill up parent, if it is not exists
        fillupParent(sline, gid, tid);

        corrRecords.push_back(libx::join(sline, '\t'));
    }
    return true;
}

void GFFCorrector::fillupParent(vector< string >& line, string& gid, string& tid)
{
    auto parent = getAttr(line.back(), PARENT);
    if (parent.empty())
    {
        if (line[2] != "mRNA" || line[2] != "exon" || line[2] != "CDS")
            return;

        auto& attrStr = line.back();
        if (!libx::endswith(attrStr, ";") && !libx::endswith(attrStr, "; "))
            attrStr.push_back(';');

        if (line[2] == "mRNA")
            attrStr += libx::join({ PARENT, gid }, kvDelim);
        else
            attrStr += libx::join({ PARENT, tid }, kvDelim);
    }
}

AnnoCorrector* createAnnoCorrector(string annoFile)
{
    if (libx::endswith(annoFile, ".gtf") || libx::endswith(annoFile, ".gtf.gz"))
        return new GTFCorrector();
    else if (libx::endswith(annoFile, ".gff") || libx::endswith(annoFile, ".gff3")
             || libx::endswith(annoFile, ".gff.gz")
             || libx::endswith(annoFile, ".gff3.gz"))
        return new GFFCorrector();
    else
        throw std::runtime_error("Invalid gtf/gff file extension: " + annoFile
                                 + ", support .gtf/.gtf.gz/.gff/.gff.gz/.gff3/.gff3.gz");
}