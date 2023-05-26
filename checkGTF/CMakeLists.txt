CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

set(SRC_PATH ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-std=c++17 -g -ggdb -W -Wall -pedantic -fopenmp -lpthread -lrt")

include_directories(../Bam2Gem)

set (src
    annoCorrector.cpp
    main.cpp
    ../Bam2Gem/geneBuilder.cpp
    ../Bam2Gem/geneFromGTF.cpp
    ../Bam2Gem/gtfReader.cpp
    ../Bam2Gem/bamUtils.cpp
    ../Bam2Gem/tagReadsWithGeneExon.cpp
)

set(EXECUTABLE_OUTPUT_PATH ${INSTALL_PATH}/bin)
add_executable(checkGTF ${src})

# CLI11 use std::filesystem, but for gcc<9,it needs a seperate library libstdc++fs
# set CLI11_HAS_FILESYSTEM=0 for prevent using std::filesystem
add_definitions(-DCLI11_HAS_FILESYSTEM=0)

target_link_libraries(checkGTF
    z
    bz2
    hts
)

install(TARGETS checkGTF RUNTIME DESTINATION bin)
set_target_properties(checkGTF PROPERTIES OUTPUT_NAME "checkGTF")
