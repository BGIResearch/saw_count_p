#include <math.h>

#include <iostream>

#include <doctest.h>

#include "utils.h"

TEST_SUITE("testing utils")
{
    TEST_CASE("Transform between (umi)sequence and unsigned long")
    {
        SUBCASE("normal")
        {
            CHECK(INT2BASE[0] == 'A');
            CHECK(INT2BASE[1] == 'C');
            CHECK(INT2BASE[2] == 'G');
            CHECK(INT2BASE[3] == 'T');
            CHECK(sizeof(INT2BASE) == (BASES_NUM + 1));

            CHECK_EQ(BASE2INT('A'), 0);
            CHECK_EQ(BASE2INT('C'), 1);
            CHECK_EQ(BASE2INT('G'), 2);
            CHECK_EQ(BASE2INT('T'), 3);
            CHECK_THROWS_AS(BASE2INT('N'), const std::runtime_error&);
        }
        SUBCASE("seq2ul")
        {
            CHECK_EQ(seq2ul("A"), 0);
            CHECK_EQ(seq2ul("TT"), 15);
            CHECK_EQ(seq2ul("ACGT"), (0 * 4 * 4 * 4 + 1 * 4 * 4 + 2 * 4 + 3));
            CHECK_EQ(seq2ul("TTTTTTTTTT"), ( unsigned long )(pow(2, 20) - 1));
            CHECK_EQ(seq2ul("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
                     ( unsigned long )(pow(2, 64) - 1));
            CHECK_THROWS(seq2ul("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTC"));
        }
        SUBCASE("ul2seq")
        {
            CHECK_EQ(ul2seq(0, 1), "A");
            CHECK_EQ(ul2seq(15, 2), "TT");
            CHECK_EQ(ul2seq(27, 4), "ACGT");
            CHECK_EQ(ul2seq(27, 3), "CGT");
            CHECK_EQ(ul2seq(27), "AAAAAAACGT");  // default length is 10
            CHECK_EQ(ul2seq(pow(2, 20) - 1, 10), "TTTTTTTTTT");
            CHECK_EQ(ul2seq(( unsigned long )(pow(2, 64) - 1), 32),
                     "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        }
    }

    TEST_CASE("Parse coorX/coorY/umi from bam record according tags")
    {
        string qname("S1000ABC");
        string contig("1");
        uint16 flag    = 0 & 16;
        int32  tid     = 0;
        int32  pos     = 100;
        uint8  mapq    = 255;
        uint64 n_cigar = 1;
        // cigar is 40M
        uint32 cigar[1] = { (0 | BAM_CMATCH) + ((40) << 4) };
        auto   mtid     = tid;
        auto   mpos     = pos;
        int32  isize    = 0;
        uint64 l_seq    = 0;
        char * seq = NULL, *qual = NULL;
        uint64 l_aux = 0;

        BamRecord bamRecord = createBamRecord();
        REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                         n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux)
                >= 0);
        updateIntTags(bamRecord, CX_TAG, 100);
        updateIntTags(bamRecord, CY_TAG, ( unsigned int )pow(2, 22));
        updateStrTags(bamRecord, UR_TAG, umi2hex(27, 4));

        unsigned int  cx, cy;
        unsigned long umi;
        parseRead(bamRecord, cx, cy, umi);
        CHECK_EQ(cx, 100);
        CHECK_EQ(cy, ( unsigned int )pow(2, 22));
        CHECK_EQ(umi, 27);

        destroyBamRecord(bamRecord);
    }

    TEST_CASE("Transform between coorX/coorY/gene and unsigned long")
    {
        unsigned long key;
        unsigned int  x, y, gene, x1, y1, gene1;
        x = 67891, y = 123456, gene = 40000;
        key = encodeKey(x, y, gene);

        REQUIRE(key == 298587905138400320);

        decodeKey(key, &x1, &y1, &gene1);
        REQUIRE_EQ(x1, x);
        REQUIRE_EQ(y1, y);
        REQUIRE_EQ(gene1, gene);

        // Test maximum values
        x = ( unsigned int )(pow(2, 22) - 1), y = x,
        gene = ( unsigned int )(pow(2, 20) - 1);
        key  = encodeKey(x, y, gene);

        REQUIRE(key == ( unsigned long )(pow(2, 64) - 1));

        decodeKey(key, &x1, &y1, &gene1);
        REQUIRE_EQ(x1, x);
        REQUIRE_EQ(y1, y);
        REQUIRE_EQ(gene1, gene);
    }

    TEST_CASE("Transform between umi sequence and hex")
    {
        string        umi("D5");
        unsigned long val = 213;

        REQUIRE_EQ(val, hex2umi(umi));
        REQUIRE_EQ(umi2hex(val, umi.size() * 2), umi);

        // Test maximum values
        umi = "FFFFFFFFFFFFFFFF";
        val = ( unsigned long )(pow(2, 64) - 1);
        REQUIRE_EQ(umi2hex(val, umi.size() * 2), umi);
        REQUIRE_EQ(val, hex2umi(umi));
    }
}