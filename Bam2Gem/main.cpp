/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */

#include "bam2gem.h"
#include "utils.h"

#include <ctime>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

#include <CLI11.hpp>
#include <libx/String.hpp>
#include <libx/Timer.hpp>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

static const std::string version = "2.3.2";

int main(int argc, char** argv)
{
    libx::Timer timer;
    // Parse the command line parameters.
    CLI::App app{ "Bam2Gem: mapping quality filter, deduplication, "
                  "set annotation, stat gene expression data." };
    app.footer("Bam2Gem version: " + version);
    app.get_formatter()->column_width(40);

    Arguments arguments;

    // Required parameters
    string bamFiles;
    app.add_option("-I,-i", bamFiles,
                   "Input bam filename or file list separated by comma")
        ->required();
    app.add_option("-O,-o", arguments.outputBamFile, "Output bam filename")->required();
    app.add_option("-A,-a", arguments.geneAnnotationFile, "Input annotation filename")
        ->check(CLI::ExistingFile)
        ->required();
    app.add_option("-S,-s", arguments.outputMetricsFile, "Output summary filename")
        ->required();
    app.add_option("-E,-e", arguments.outputExpFile,
                   "Output barcode gene expression filename")
        ->required();
    app.add_option("--sn", arguments.sn, "STOmics Chip Serial Number")->required();

    // Optional parameters
    arguments.mapQualThre = 10;
    app.add_option("-Q,-q", arguments.mapQualThre,
                   "Set mapping quality threshold, default 10")
        ->check(CLI::PositiveNumber);
    arguments.coreNum = std::thread::hardware_concurrency();
    app.add_option("-C,-c", arguments.coreNum, "Set cpu cores, default detect")
        ->check(CLI::PositiveNumber);
    arguments.memoryGB = getSystemMemoryGB();
    app.add_option("-M,-m", arguments.memoryGB,
                   "Set avaliable memory(GB), default detect")
        ->check(CLI::PositiveNumber);
    arguments.saveLowQualReads = false;
    app.add_flag("--save_lq", arguments.saveLowQualReads,
                 "Save low quality reads, default false");
    arguments.saveDupReads = false;
    app.add_flag("--save_dup", arguments.saveDupReads,
                 "Save duplicate reads, default false");
    // Umi correction parameters
    arguments.umiPara.enable = false;
    app.add_flag("--umi_on", arguments.umiPara.enable,
                 "Enable umi correction, default disable");
    arguments.umiPara.minUmiNum = 5;
    app.add_option("--umi_min_num", arguments.umiPara.minUmiNum,
                   "Minimum umi number for correction, default 5")
        ->check(CLI::PositiveNumber);
    arguments.umiPara.mismatch = 1;
    app.add_option("--umi_mismatch", arguments.umiPara.mismatch,
                   "Maximum mismatch for umi correction, default 1")
        ->check(CLI::PositiveNumber);
    arguments.umiPara.len = 10;
    app.add_option("--umi_len", arguments.umiPara.len, "UMI length, default 10")
        ->check(CLI::PositiveNumber);
    app.add_option("--sat_file", arguments.outputSatFile,
                   "Output sequencing saturation file, default None");
    arguments.multiMap = false;
    app.add_flag("--multi_map", arguments.multiMap,
                 "Enable multi-mapping reads correction, default disable");

    CLI11_PARSE(app, argc, argv);

    // Check the input bam files
    vector< string > bamList;
    libx::split(bamFiles, ',', bamList);
    if (bamList.empty())
    {
        std::cerr << ERR_CODE << "001 "
                  << "Invalid parameter of -I: " << bamFiles << std::endl;
        exit(-1);
    }
    else
    {
        for (auto& f : bamList)
        {
            if (!fs::exists(f))
            {
                std::cerr << ERR_CODE << "001 "
                          << "Not exists bam file: " << f << std::endl;
                exit(-1);
            }
        }
    }
    arguments.inputBamFiles = bamList;

    // Check output parameters
    for (auto f : { arguments.outputBamFile, arguments.outputExpFile,
                    arguments.outputMetricsFile, arguments.outputSatFile })
    {
        if (f.empty())
            continue;
        if (!fs::exists(fs::absolute(fs::path(f)).parent_path()))
        {
            std::cerr << ERR_CODE << "001 "
                      << "Not exists parent path of specific file: " << f << std::endl;
            exit(-1);
        }
    }

    // Set the default logger to file logger.
    std::time_t t =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream ostr;
    ostr << "Bam2Gem_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".log";
    try
    {
        // auto file_logger = spdlog::basic_logger_mt("main", "logs/" +
        // ostr.str());
        auto file_sink =
            std::make_shared< spdlog::sinks::basic_file_sink_mt >("logs/" + ostr.str());
        auto main_logger    = std::make_shared< spdlog::logger >("main", file_sink);
        auto gtf_logger     = std::make_shared< spdlog::logger >("gtf", file_sink);
        auto process_logger = std::make_shared< spdlog::logger >("process", file_sink);
        spdlog::register_logger(main_logger);
        spdlog::register_logger(gtf_logger);
        spdlog::register_logger(process_logger);
        spdlog::set_default_logger(process_logger);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cerr << ERR_CODE << "001 "
                  << "Log init failed: " << ex.what() << std::endl;
        exit(-2);
    }
    spdlog::set_level(spdlog::level::debug);  // Set global log level.
    spdlog::flush_on(spdlog::level::debug);
    spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %L %n: %v");

    spdlog::get("main")->info("{} {}", argv[0], arguments.str());

    std::exception_ptr ex;
    try
    {
        Bam2Gem bam2gem(arguments);
        bam2gem.prepare();
        bam2gem.doWork();
    }
    catch (const std::exception& e)
    {
        std::cerr << ERR_CODE << "003 " << e.what() << std::endl;
        ex = std::current_exception();
        std::rethrow_exception(ex);
    }

    spdlog::get("main")->info("Bam2Gem done. Elapsed time(s):{}", timer.toc());

    return 0;
}
