/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#pragma once

#include <map>
#include <string>
#include <vector>
using namespace std;

class AnnoCorrector
{
public:
    explicit AnnoCorrector();
    virtual ~AnnoCorrector();

    bool         load(string filename);
    virtual bool correction() = 0;
    bool         dump(string filename);

protected:
    virtual void   checkStrand(string& strand);
    virtual string fillupAttr(string& attr);
    virtual string getAttr(string& attr, string key);

protected:
    vector< string >      rawRecords, corrRecords;
    char                  kvDelim;  // Delimiter of key-value pairs in attribute column
    map< string, string > keyMap;
    string                gidName;
};

class GTFCorrector : public AnnoCorrector
{
public:
    explicit GTFCorrector();

    virtual bool correction();

protected:
    vector< string > fillupLine(vector< vector< string > >& cache);

private:
    inline static const string GID   = "gene_id";
    inline static const string GNAME = "gene_name";
    inline static const string TID   = "transcript_id";
    inline static const string TNAME = "transcript_name";
};

class GFFCorrector : public AnnoCorrector
{
public:
    explicit GFFCorrector();

    virtual bool correction();

protected:
    void fillupParent(vector< string >& line, string& gid, string& tid);

private:
    inline static const string ID     = "ID";
    inline static const string NAME   = "Name";
    inline static const string PARENT = "Parent";
};

AnnoCorrector* createAnnoCorrector(string annoFile);