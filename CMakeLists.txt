cmake_minimum_required(VERSION 3.5)


set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TURE)

project(RNASeq)

#set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

if ( UNITTEST STREQUAL "OFF" )
    # the main project
    add_subdirectory(Bam2Gem)
    add_subdirectory(checkGTF)
    message(STATUS "Compile Bam2Gem")
else()
    # the unittest project
    add_subdirectory(unittest)
    message(STATUS "Compile unittest")
endif()

