
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

set(INCLUDE_PATH ${PROJECT_SOURCE_DIR}/Bam2Gem)
set(SOURCE_PATH ${PROJECT_SOURCE_DIR}/unittest)

# set(CMAKE_CXX_FLAGS "-std=c++17 -g -O0 --coverage")
# set(CMAKE_CXX_FLAGS "-std=c++17 -O2")
set(CMAKE_CXX_FLAGS "-std=c++17 -g -O0")


include_directories(${INCLUDE_PATH})

aux_source_directory(. src)
# set(src DIR_SRCS)

set(EXECUTABLE_OUTPUT_PATH ${INSTALL_PATH}/bin)

add_executable(app_tests ${src}
    ../Bam2Gem/gtfReader.cpp
    ../Bam2Gem/geneBuilder.cpp
    ../Bam2Gem/geneFromGTF.cpp
    ../Bam2Gem/tagReadsWithGeneExon.cpp
    ../Bam2Gem/bamUtils.cpp
    ../Bam2Gem/utils.cpp
)

target_link_libraries(app_tests
    z
    bz2
    hts
)

set_target_properties(app_tests PROPERTIES OUTPUT_NAME "unittest")