#include <doctest.h>
#include <htslib/sam.h>
#include <libx/File.hpp>
#include <libx/GZFile.hpp>
#include <libx/System.hpp>

#include <string>
using namespace std;

#include "annotationException.h"
#include "bamUtils.h"
#include "geneBuilder.h"
#include "gtfReader.h"
#include "tagReadsWithGeneExon.h"
#include "utils.h"

// clang-format off
// gtf format with one gene, one transcript and two exons
constexpr char gtfRecords[] =
    "1\tensembl_havana\tgene\t3205901\t3671498\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; gene_version \"5\"; gene_name \"Xkr4\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"\n"
    "1\thavana\ttranscript\t3205901\t3216344\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; gene_version \"5\"; transcript_id \"ENSMUST00000162897\"; transcript_version \"1\"; gene_name \"Xkr4\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"Xkr4-203\"; transcript_source \"havana\"; transcript_biotype \"processed_transcript\"; transcript_support_level \"1\"\n"
    "1\thavana\texon\t3213609\t3216344\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; gene_version \"5\"; transcript_id \"ENSMUST00000162897\"; transcript_version \"1\"; exon_number \"1\"; gene_name \"Xkr4\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"Xkr4-203\"; transcript_source \"havana\"; transcript_biotype \"processed_transcript\"; exon_id \"ENSMUSE00000858910\"; exon_version \"1\"; transcript_support_level \"1\"\n"
    "1\thavana\texon\t3205901\t3207317\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; gene_version \"5\"; transcript_id \"ENSMUST00000162897\"; transcript_version \"1\"; exon_number \"2\"; gene_name \"Xkr4\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"Xkr4-203\"; transcript_source \"havana\"; transcript_biotype \"processed_transcript\"; exon_id \"ENSMUSE00000866652\"; exon_version \"1\"; transcript_support_level \"1\"\n";
// gff format with one gene, one transcript and three exons
constexpr char gffRecords[] =
    "##gff-version\t3\nChr1\tAraport11\tgene\t3631\t5899\t.\t+\t.\tID=AT1G01010;Name=AT1G01010\n"
    "Chr1\tAraport11\tmRNA\t3631\t5899\t.\t+\t.\tID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1\n"
    "Chr1\tAraport11\tfive_prime_UTR\t3631\t3759\t.\t+\t.\tID=AT1G01010:five_prime_UTR:1;Parent=AT1G01010.1;Name=NAC001:five_prime_UTR:1\n"
    "Chr1\tAraport11\texon\t3631\t3913\t.\t+\t.\tID=AT1G01010:exon:1;Parent=AT1G01010.1;Name=NAC001:exon:1\n"
    "Chr1\tAraport11\tCDS\t3760\t3913\t.\t+\t0\tID=AT1G01010:CDS:1;Parent=AT1G01010.1;Name=NAC001:CDS:1\n"
    "Chr1\tAraport11\tprotein\t3760\t5630\t.\t+\t.\tID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1\n"
    "Chr1\tAraport11\texon\t3996\t4276\t.\t+\t.\tID=AT1G01010:exon:2;Parent=AT1G01010.1;Name=NAC001:exon:2\n"
    "Chr1\tAraport11\tCDS\t3996\t4276\t.\t+\t2\tID=AT1G01010:CDS:2;Parent=AT1G01010.1;Name=NAC001:CDS:2\n"
    "Chr1\tAraport11\texon\t4486\t4605\t.\t+\t.\tID=AT1G01010:exon:3;Parent=AT1G01010.1;Name=NAC001:exon:3\n";
// gtf format with three genes, first gene has one transcript and two exons, second gene has two transcripts with one and two exons
// the third gene is copied from gene1 but has negative strand
/*
 * gene1  100                       500
 * tran1  100                   450
 * exon1  100     200
 * exon2                    400 450
 *
 * gene2    150                                     800
 * tran1    150              420
 * exon1    150              420
 * tran2      180                            600
 * exon1      180      250
 * exon2                  300       500
 * exon3                                550  600
 *
 * gene3  100                       500
 * tran1  100                   450
 * exon1  100     200
 * exon2                    400 450
*/
constexpr char multiGeneGtfRecords[] = 
    "1\tensembl_havana\tgene\t100\t500\t.\t+\t.\tgene_id \"gid1\"; gene_name \"gene1\"\n"
    "1\thavana\ttranscript\t100\t450\t.\t+\t.\tgene_id \"gid1\"; transcript_id \"tid1\"; gene_name \"gene1\";  transcript_name \"trans1\";\n"
    "1\thavana\texon\t100\t200\t.\t+\t.\tgene_id \"gid1\"; transcript_id \"tid1\"; gene_name \"gene1\"; transcript_name \"trans1\"; exon_id \"eid1\"\n"
    "1\thavana\texon\t400\t450\t.\t+\t.\tgene_id \"gid1\"; transcript_id \"tid1\"; gene_name \"gene1\"; transcript_name \"trans1\"; exon_id \"eid2\"\n"
    "1\tensembl_havana\tgene\t150\t800\t.\t-\t.\tgene_id \"gid2\"; gene_name \"gene2\"\n"
    "1\thavana\ttranscript\t150\t420\t.\t-\t.\tgene_id \"gid2\"; transcript_id \"tid1\"; gene_name \"gene2\";  transcript_name \"trans1\";\n"
    "1\thavana\texon\t150\t420\t.\t-\t.\tgene_id \"gid2\"; transcript_id \"tid1\"; gene_name \"gene2\"; transcript_name \"trans1\"; exon_id \"eid1\"\n"
    "1\thavana\ttranscript\t180\t600\t.\t-\t.\tgene_id \"gid2\"; transcript_id \"tid2\"; gene_name \"gene2\";  transcript_name \"trans2\";\n"
    "1\thavana\texon\t180\t250\t.\t-\t.\tgene_id \"gid2\"; transcript_id \"tid2\"; gene_name \"gene2\"; transcript_name \"trans2\"; exon_id \"eid1\"\n"
    "1\thavana\texon\t300\t500\t.\t-\t.\tgene_id \"gid2\"; transcript_id \"tid2\"; gene_name \"gene2\"; transcript_name \"trans2\"; exon_id \"eid2\"\n"
    "1\thavana\texon\t550\t600\t.\t-\t.\tgene_id \"gid2\"; transcript_id \"tid2\"; gene_name \"gene2\"; transcript_name \"trans2\"; exon_id \"eid3\"\n"
    "1\tensembl_havana\tgene\t100\t500\t.\t-\t.\tgene_id \"gid3\"; gene_name \"gene3\"\n"
    "1\thavana\ttranscript\t100\t450\t.\t-\t.\tgene_id \"gid3\"; transcript_id \"tid1\"; gene_name \"gene3\";  transcript_name \"trans1\";\n"
    "1\thavana\texon\t100\t200\t.\t-\t.\tgene_id \"gid3\"; transcript_id \"tid1\"; gene_name \"gene3\"; transcript_name \"trans1\"; exon_id \"eid1\"\n"
    "1\thavana\texon\t400\t450\t.\t-\t.\tgene_id \"gid3\"; transcript_id \"tid1\"; gene_name \"gene3\"; transcript_name \"trans1\"; exon_id \"eid2\"\n";
// clang-format on

TEST_SUITE("testing gtfReader")
{
    TEST_CASE("normal gtf/gff format")
    {
        AnnoReader* annoReader = nullptr;

        SUBCASE("gtf")
        {
            for (int i = 0; i < 2; ++i)
            {
                string filename;
                if (i == 0)
                {
                    filename = "normalTest.gtf";
                    libx::writeFile(gtfRecords, filename);
                }
                else
                {
                    filename = "normalTest.gtf.gz";
                    libx::writeGZFile(gtfRecords, filename);
                }

                unordered_map< string, vector< GTFRecord > > gatherByGeneName;
                annoReader = createAnnoReader(filename);
                annoReader->loadAnnoFile(gatherByGeneName);

                CHECK(gatherByGeneName.size() == 1);
                REQUIRE(gatherByGeneName.count("Xkr4"));
                auto& vec = gatherByGeneName["Xkr4"];
                CHECK(vec.size() == 4);
                CHECK(vec[0].getGeneName() == "Xkr4");
                CHECK(vec[0].getGeneID() == "ENSMUSG00000051951");
                CHECK(vec[1].getGeneName() == "Xkr4");
                CHECK(vec[1].getGeneID() == "ENSMUSG00000051951");
                CHECK(vec[1].getTranscriptName() == "Xkr4-203");
                CHECK(vec[1].getTranscriptID() == "ENSMUST00000162897");
                CHECK(vec[0].getFeatureType() == "gene");
                CHECK(vec[0].getStart() == 3205901);
                CHECK(vec[0].getEnd() == 3671498);
                CHECK(vec[0].getContig() == "1");
                CHECK(vec[0].isNegativeStrand());

                if (annoReader != nullptr)
                {
                    delete annoReader;
                    annoReader = nullptr;
                }
                // Remove temp file
                string removeCmd("rm " + filename);
                REQUIRE(libx::subprocess(removeCmd) == 0);
            }
        }
        SUBCASE("gff")
        {
            for (int i = 0; i < 4; ++i)
            {
                string filename;
                if (i < 2)
                    filename = "normalTest.gff";
                else
                    filename = "normalTest.gff3";

                if (i % 2 == 0)
                {
                    libx::writeFile(gffRecords, filename);
                }
                else
                {
                    filename += ".gz";
                    libx::writeGZFile(gffRecords, filename);
                }

                unordered_map< string, vector< GTFRecord > > gatherByGeneName;
                annoReader = createAnnoReader(filename);
                annoReader->loadAnnoFile(gatherByGeneName);

                CHECK(gatherByGeneName.size() == 1);
                REQUIRE(gatherByGeneName.count("AT1G01010"));
                auto& vec = gatherByGeneName["AT1G01010"];
                CHECK(vec.size() == 8);
                CHECK(vec[0].getGeneName() == "AT1G01010");
                CHECK(vec[0].getGeneID() == "AT1G01010");
                CHECK(vec[7].getGeneName() == "AT1G01010");
                CHECK(vec[7].getGeneID() == "AT1G01010");
                CHECK(vec[7].getTranscriptName() == "AT1G01010.1");
                CHECK(vec[7].getTranscriptID() == "AT1G01010.1");
                CHECK(vec[0].getFeatureType() == "gene");
                CHECK(vec[0].getStart() == 3631);
                CHECK(vec[0].getEnd() == 5899);
                CHECK(vec[0].getContig() == "Chr1");
                CHECK(!vec[0].isNegativeStrand());

                if (annoReader != nullptr)
                {
                    delete annoReader;
                    annoReader = nullptr;
                }
                // Remove temp file
                string removeCmd("rm " + filename);
                REQUIRE(libx::subprocess(removeCmd) == 0);
            }
        }
    }

    TEST_CASE("invalid file format")
    {
        AnnoReader*                                  annoReader = nullptr;
        unordered_map< string, vector< GTFRecord > > gatherByGeneName;
        string                                       filename;

        SUBCASE("invalid file name suffix")
        {
            REQUIRE_THROWS_AS(annoReader = createAnnoReader("test.gff4"),
                              const AnnotationException&);
            REQUIRE_THROWS_AS(annoReader = createAnnoReader("test"),
                              const AnnotationException&);
            REQUIRE_THROWS_AS(annoReader = createAnnoReader("test.gff.gtf3"),
                              const AnnotationException&);
            REQUIRE_THROWS_AS(annoReader = createAnnoReader("test.gz"),
                              const AnnotationException&);

            CHECK_NOTHROW(annoReader = createAnnoReader("test.gtf"));
            CHECK_NOTHROW(annoReader = createAnnoReader("test.gtf.gz"));
            CHECK_NOTHROW(annoReader = createAnnoReader("test.gff"));
            CHECK_NOTHROW(annoReader = createAnnoReader("test.gff.gz"));
            CHECK_NOTHROW(annoReader = createAnnoReader("test.gff3"));
            CHECK_NOTHROW(annoReader = createAnnoReader("test.gff3.gz"));
        }
        SUBCASE("invalid range")
        {
            // 2^31 = 2147483648
            string records =
                "1\tensembl_havana\tgene\t3205901\t2147483648\t.\t-\t.\tgene_"
                "id \"ENSMUSG00000051951\"; gene_name \"Xkr4\"\n";
            filename = "invalid_range.gtf";
            libx::writeFile(records, filename);

            unordered_map< string, vector< GTFRecord > > gatherByGeneName;
            annoReader = createAnnoReader(filename);
            REQUIRE_THROWS_AS(annoReader->loadAnnoFile(gatherByGeneName), exception&);
        }
        SUBCASE("lack of attributes in gtf")
        {
            // four lines are invalid, one line is valid
            // clang-format off
            string records =
                "1\tensembl_havana\tgene\t3205901\t3205903\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"\n"
                "1\tensembl_havana\tgene\t3205901\t3205903\t.\t-\t.\tgene_name \"Xkr4\"\n"
                "1\thavana\ttranscript\t3205901\t3216344\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; transcript_id \"ENSMUST00000162897\"; gene_name \"Xkr4\"\n"
                "1\thavana\ttranscript\t3205901\t3216344\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; gene_name \"Xkr4\"; transcript_name \"Xkr4-203\"\n"
                "1\thavana\ttranscript\t3205901\t3216344\t.\t-\t.\tgene_id \"ENSMUSG00000051951\"; transcript_id \"ENSMUST00000162897\"; gene_name \"Xkr4\"; transcript_name \"Xkr4-203\"\n";
            // clang-format on
            filename = "lack_attributes.gtf";
            libx::writeFile(records, filename);

            unordered_map< string, vector< GTFRecord > > gatherByGeneName;
            annoReader = createAnnoReader(filename);
            annoReader->loadAnnoFile(gatherByGeneName);
            CHECK(gatherByGeneName.size() == 1);
            CHECK(gatherByGeneName["Xkr4"].size() == 1);
        }
        SUBCASE("lack of attributes in gff")
        {
            // five lines are invalid, two line is valid
            // clang-format off
            string records =
                "Chr1\tAraport11\tgene\t3631\t5899\t.\t+\t.\tID=AT1G01010\n"
                "Chr1\tAraport11\tgene\t3631\t5899\t.\t+\t.\tName=AT1G01010\n"
                "Chr1\tAraport11\tmRNA\t3631\t5899\t.\t+\t.\tParent=AT1G01010;Name=AT1G01010.1\n"
                "Chr1\tAraport11\tmRNA\t3631\t5899\t.\t+\t.\tID=AT1G01010.1;Parent=AT1G01010\n"
                "Chr1\tAraport11\tmRNA\t3631\t5899\t.\t+\t.\tID=AT1G01010.1;Name=AT1G01010.1\n"
                "Chr1\tAraport11\tgene\t3631\t5899\t.\t+\t.\tID=AT1G01010;Name=AT1G01010\n"
                "Chr1\tAraport11\tmRNA\t3631\t5899\t.\t+\t.\tID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1\n";
            // clang-format on
            filename = "lack_attributes.gff";
            libx::writeFile(records, filename);

            unordered_map< string, vector< GTFRecord > > gatherByGeneName;
            annoReader = createAnnoReader(filename);
            annoReader->loadAnnoFile(gatherByGeneName);
            CHECK(gatherByGeneName.size() == 1);
            CHECK(gatherByGeneName["AT1G01010"].size() == 2);
        }

        if (annoReader != nullptr)
        {
            delete annoReader;
            annoReader = nullptr;
        }
        // Remove temp file
        if (!filename.empty())
        {
            string removeCmd("rm " + filename);
            REQUIRE(libx::subprocess(removeCmd) == 0);
        }
    }
}

TEST_SUITE("testing GeneBuilder")
{
    TEST_CASE("normal GeneBuilder")
    {
        AnnoReader* annoReader = nullptr;

        string filename("geneBuilderTest.gtf");
        libx::writeFile(gtfRecords, filename);

        unordered_map< string, vector< GTFRecord > > gatherByGeneName;
        annoReader = createAnnoReader(filename);
        annoReader->loadAnnoFile(gatherByGeneName);

        GeneBuilder           geneBuilder;
        GTFIterator           gtfIterator = gatherByGeneName.begin();
        vector< GeneFromGTF > genes       = geneBuilder.makeGene(gtfIterator);
        REQUIRE(genes.size() == 1);
        auto& gene = genes.front();
        REQUIRE(gene.isNegativeStrand());
        REQUIRE(gene.getContig() == "1");
        REQUIRE(gene.getName() == "Xkr4");
        REQUIRE(gene.getStart() == 3205901);
        REQUIRE(gene.getEnd() == 3671498);
        auto& transcriptMap = gene.getTranscripts();
        REQUIRE(transcriptMap.size() == 1);
        REQUIRE(transcriptMap.count("ENSMUST00000162897") == 1);
        auto& transcript = transcriptMap.at("ENSMUST00000162897");
        REQUIRE(transcript.getTranscriptionStart() == 3205901);
        REQUIRE(transcript.getTranscriptionEnd() == 3216344);
        auto& exons = transcript.getExons();
        REQUIRE(exons.size() == 2);
        // Exons stored by start order
        int s1 = 3205901, e1 = 3207317;
        int s2 = 3213609, e2 = 3216344;
        REQUIRE(exons[0].start == s1);
        REQUIRE(exons[0].end == e1);
        REQUIRE(exons[1].start == s2);
        REQUIRE(exons[1].end == e2);
        // Testing method of class TranscriptFromGTF
        CHECK(!transcript.inExon(s1 - 1));
        CHECK(transcript.inExon(s1));
        CHECK(transcript.inExon(s1 + 1));
        CHECK(transcript.inExon(e2 - 1));
        CHECK(transcript.inExon(e2));
        CHECK(!transcript.inExon(e2 + 1));
        CHECK(transcript.inExon(e1));
        CHECK(!transcript.inExon(e1 + 1));

        pair< int, int > max_cnts{ 0, 0 };
        transcript.assignLocusFunction(s1, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(100, 0));
        max_cnts = { 0, 0 };
        transcript.assignLocusFunction(s1 - 50, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(50, 0));
        transcript.assignLocusFunction(e1 - 49, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(50, 50));
        transcript.assignLocusFunction(e1 - 89, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(90, 10));
        transcript.assignLocusFunction(e1 - 89, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(90, 10));
        transcript.assignLocusFunction(s2 - 8, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(92, 8));
        transcript.assignLocusFunction(e2 - 98, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(99, 0));
        transcript.assignLocusFunction(s2 - 1, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(99, 1));
        transcript.assignLocusFunction(e2 - 99, 100, max_cnts);
        REQUIRE(max_cnts == make_pair(100, 0));

        if (annoReader != nullptr)
        {
            delete annoReader;
            annoReader = nullptr;
        }
        // Remove temp file
        string removeCmd("rm " + filename);
        REQUIRE(libx::subprocess(removeCmd) == 0);
    }

    TEST_CASE("multi genes")
    {
        AnnoReader* annoReader = nullptr;

        string filename("geneBuilderTest.gtf");
        libx::writeFile(multiGeneGtfRecords, filename);

        unordered_map< string, vector< GTFRecord > > gatherByGeneName;
        annoReader = createAnnoReader(filename);
        annoReader->loadAnnoFile(gatherByGeneName);

        GeneBuilder geneBuilder;
        GTFIterator gtfIterator = gatherByGeneName.begin();
        for (; gtfIterator != gatherByGeneName.end(); ++gtfIterator)
        {
            vector< GeneFromGTF > genes = geneBuilder.makeGene(gtfIterator);
            if (gtfIterator->first == "gene1")
            {
                REQUIRE(genes.size() == 1);
                auto& gene = genes.front();
                REQUIRE(gene.getName() == gtfIterator->first);
                auto& transcriptMap = gene.getTranscripts();
                REQUIRE(transcriptMap.size() == 1);
                auto& transcript = transcriptMap.at("tid1");
                CHECK(transcript.getTranscriptName() == "trans1");
                auto& exons = transcript.getExons();
                REQUIRE(exons.size() == 2);
            }
            else if (gtfIterator->first == "gene2")
            {
                REQUIRE(genes.size() == 1);
                auto& gene = genes.front();
                REQUIRE(gene.getName() == gtfIterator->first);
                auto& transcriptMap = gene.getTranscripts();
                REQUIRE(transcriptMap.size() == 2);
                auto transcript = transcriptMap.at("tid1");
                CHECK(transcript.getTranscriptName() == "trans1");
                auto exons = transcript.getExons();
                REQUIRE(exons.size() == 1);
                transcript = transcriptMap.at("tid2");
                CHECK(transcript.getTranscriptName() == "trans2");
                exons = transcript.getExons();
                REQUIRE(exons.size() == 3);
            }
        }

        if (annoReader != nullptr)
        {
            delete annoReader;
            annoReader = nullptr;
        }
        // Remove temp file
        string removeCmd("rm " + filename);
        REQUIRE(libx::subprocess(removeCmd) == 0);
    }
}

TEST_SUITE("testing TagReadsWithGeneExon")
{
    TEST_CASE("normal")
    {
        string filename("TagReadsWithGeneExonTest.gtf");
        libx::writeFile(multiGeneGtfRecords, filename);

        TagReadsWithGeneExon tagReads(filename);
        tagReads.makeOverlapDetector();

        BamRecord bamRecord = createBamRecord();
        string    qname("S1000ABC");
        string    contig("1");
        uint16    flag    = 0 | 16;
        int32     tid     = 0;
        int32     pos     = 100;
        uint8     mapq    = 255;
        uint64    n_cigar = 1;
        // cigar is 40M
        uint32 cigar[1] = { (0 | BAM_CMATCH) + ((40) << 4) };
        auto   mtid     = tid;
        auto   mpos     = pos;
        int32  isize    = 0;
        uint64 l_seq    = 0;
        char * seq = NULL, *qual = NULL;
        uint64 l_aux = 0;

        string       outGene;
        unsigned int outStrand = 0;
        unsigned int outLocus  = 0;
        // start 100 cigar 40M
        SUBCASE("mapped to one gene and one transcript, total mapped to exon")
        {
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                             n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene3");
            CHECK(outStrand == 2);
            CHECK(outLocus == 3);
        }
        // start 79 cigar 40M
        SUBCASE("mapped to one gene and one transcript, half mapped to exon")
        {
            pos = 79;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                             n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene3");
            CHECK(outStrand == 2);
            CHECK(outLocus == 3);
        }
        // start 78 cigar 40M, convert 0-based to 1-based
        SUBCASE("mapped to one gene and one transcript, less mapped to exon")
        {
            pos = 78;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                             n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "");
            CHECK(outStrand == 0);
            CHECK(outLocus == 1);
        }
        // start 20 cigar 40M
        SUBCASE("no mapped")
        {
            pos = 20;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                             n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "");
            CHECK(outStrand == 0);
            CHECK(outLocus == 1);
        }
        // start 501 cigar 40M
        SUBCASE("mapped to one gene and one transcript, total mapped to intron")
        {
            pos = 501;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 2);
        }
        // start 531 cigar 40M
        SUBCASE("mapped to one gene and one transcript, exon > intron")
        {
            pos = 531;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 3);
        }
        // start 530 cigar 40M
        SUBCASE("mapped to one gene and one transcript, exon == intron")
        {
            pos = 530;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 3);
        }
        // start 529 cigar 40M
        SUBCASE("mapped to one gene and one transcript, exon < intron")
        {
            pos = 529;
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 2);
        }

        // start 400 cigar 80M
        SUBCASE("mapped to two genes, exon > intron")
        {
            pos      = 400;
            cigar[0] = { (0 | BAM_CMATCH) + ((80) << 4) };
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 3);
        }
        // start 451 cigar 20M30N20M
        SUBCASE("multi align blocks, exon")
        {
            pos              = 451;
            n_cigar          = 3;
            uint32 cigar2[3] = { (0 | BAM_CMATCH) + ((20) << 4),
                                 (0 | BAM_CREF_SKIP) + ((30) << 4),
                                 (0 | BAM_CMATCH) + ((20) << 4) };
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar2, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 3);
        }
        // start 451 cigar 19M31N20M
        SUBCASE("multi align blocks, intron")
        {
            pos              = 451;
            n_cigar          = 3;
            uint32 cigar2[3] = { (0 | BAM_CMATCH) + ((19) << 4),
                                 (0 | BAM_CREF_SKIP) + ((31) << 4),
                                 (0 | BAM_CMATCH) + ((20) << 4) };
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar2, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "gene2");
            CHECK(outStrand == 2);
            CHECK(outLocus == 2);
        }
        // start 451 cigar 39M10000N40M
        SUBCASE("multi align blocks, intergenic")
        {
            pos              = 451;
            n_cigar          = 3;
            uint32 cigar2[3] = { (0 | BAM_CMATCH) + ((39) << 4),
                                 (0 | BAM_CREF_SKIP) + ((10000) << 4),
                                 (0 | BAM_CMATCH) + ((40) << 4) };
            REQUIRE(bam_set1(bamRecord, qname.size(), qname.c_str(), flag, tid, pos - 1,
                             mapq, n_cigar, cigar2, mtid, mpos, isize, l_seq, seq, qual,
                             l_aux)
                    >= 0);

            tagReads.setAnnotation(bamRecord, contig, outGene, &outStrand, &outLocus);
            CHECK(outGene == "");
            CHECK(outStrand == 0);
            CHECK(outLocus == 1);
        }

        destroyBamRecord(bamRecord);
        // Remove temp file
        string removeCmd("rm " + filename);
        REQUIRE(libx::subprocess(removeCmd) == 0);
    }
}