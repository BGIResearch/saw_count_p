#include <doctest.h>

#include "umiCorrection.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

TEST_SUITE("testing umi correction")
{
    TEST_CASE("normal umi correction")
    {
        typedef unsigned long                                    mytype;
        constexpr unsigned int                                   umilen = 3;
        UmiCorrection< mytype >                                  umiCorrection;
        unordered_map< mytype, unordered_map< mytype, int > >    umiMismatch;
        unordered_map< mytype, unordered_map< mytype, mytype > > umiCorrect;

        // prepare umis and counts
        vector< string > seqs = { "AAA", "GGA", "AGA", "AAT", "GGG", "CCC" };
        vector< mytype > umis(seqs.size());
        std::transform(seqs.begin(), seqs.end(), umis.begin(),
                       [&](string& s) { return seq2ul(s); });
        vector< int > cnts = { 5, 4, 3, 2, 1, 1 };
        umiMismatch[1]     = {};
        for (int i = 0; i < umis.size(); ++i)
        {
            mytype umi          = umis[i];
            auto   cnt          = cnts[i];
            umiMismatch[1][umi] = cnt;
        }

        // deduplication of umis
        umiCorrection.deDupUmi(umiMismatch, umiCorrect);

        // check the results
        unordered_map< string, string > dict;
        for (auto& [_, pairs] : umiCorrect)
        {
            for (auto& [u1, u2] : pairs)
                dict[ul2seq(u1, umilen)] = ul2seq(u2, umilen);
        }
        CHECK(dict.size() == 3);
        CHECK(dict["GGG"] == "GGA");
        CHECK(dict["AGA"] == "AAA");
        CHECK(dict["AAT"] == "AAA");

        auto& res = umiMismatch[1];
        CHECK(res[umis[0]] == (5 + 3 + 2));
        CHECK(res[umis[1]] == (4 + 1));
        CHECK(res[umis[2]] == (0));
        CHECK(res[umis[3]] == (0));
        CHECK(res[umis[4]] == (0));
        CHECK(res[umis[5]] == (1));
    }
}