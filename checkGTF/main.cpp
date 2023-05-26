/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include "annoCorrector.h"
#include "tagReadsWithGeneExon.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
namespace fs = std::filesystem;
using namespace std;

#include <CLI11.hpp>
#include <libx/String.hpp>
#include <libx/Timer.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

int check(string gtfFile)
{
    // Load annotations.
    TagReadsWithGeneExon tagReadsWithGeneExon(gtfFile);
    if (tagReadsWithGeneExon.makeOverlapDetector() != 0)
    {
        throw std::runtime_error(
            "Failed makeOverlapDetector! Please check annotation file format.");
    }
    // The first string is empty
    vector< string > genes = tagReadsWithGeneExon.getGeneName();
    return genes.size() - 1;
}

void update(string inputFile, string outputFile)
{
    AnnoCorrector* corr = createAnnoCorrector(inputFile);
    corr->load(inputFile);
    corr->correction();
    corr->dump(outputFile);

    if (corr != nullptr)
    {
        delete corr;
        corr = nullptr;
    }
}

int main(int argc, char** argv)
{
    libx::Timer timer;
    // Parse the command line parameters.
    CLI::App app;

    // Required parameters
    string inputGtf, outputGtf;
    app.add_option("-I,-i", inputGtf, "Input gtf/gff filename")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-O,-o", outputGtf, "Output gtf/gff filename after formatting");

    CLI11_PARSE(app, argc, argv);

    // Set the default logger to file logger.
    try
    {
        auto console_sink   = std::make_shared< spdlog::sinks::stdout_color_sink_mt >();
        auto process_logger = std::make_shared< spdlog::logger >("process", console_sink);
        auto gtf_logger     = std::make_shared< spdlog::logger >("gtf", console_sink);
        spdlog::register_logger(process_logger);
        spdlog::register_logger(gtf_logger);
        spdlog::set_default_logger(process_logger);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cerr << "Log init failed: " << ex.what() << std::endl;
        exit(-1);
    }
    spdlog::set_level(spdlog::level::debug);  // Set global log level.
    spdlog::flush_on(spdlog::level::debug);
    spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %L %n: %v");

    int rawGenes = check(inputGtf);
    if (!outputGtf.empty())
    {
        update(inputGtf, outputGtf);
        int newGenes = check(outputGtf);
        spdlog::info("Rescue {} genes", newGenes - rawGenes);
    }

    spdlog::info("checkGTF done. Elapsed time(s):{}", timer.toc());

    return 0;
}
