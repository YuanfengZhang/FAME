cmake_minimum_required(VERSION 3.12)
project(FAME VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags for optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find dependencies
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

# Main executable
add_executable(fame
    main.cpp
    RefGenome.cpp
    ReadQueue.cpp
    Read.cpp
    DnaBitStr.cpp
    ShiftAnd.cpp
    LevenshtDP.cpp
    CONST.cpp
    RefReader.cpp
    RefReader_istr.cpp
    gzstream/gzstream.C
)

target_include_directories(fame PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/sparsehash/include/usr/local/include
    ${CMAKE_SOURCE_DIR}/hopscotch-map/include/tsl
    ${CMAKE_SOURCE_DIR}/spaced_nthash
)

target_link_libraries(fame PRIVATE
    OpenMP::OpenMP_CXX
    ZLIB::ZLIB
)

# Installation
install(TARGETS fame DESTINATION bin)

# Optional: install headers for development
install(FILES
    CONST.h
    ReadQueue.h
    RefGenome.h
    Read.h
    DnaBitStr.h
    ShiftAnd.h
    LevenshtDP.h
    RefReader.h
    RefReader_istr.h
    structs.h
    BandedMatrix.h
    BitFunctions.h
    DESTINATION include/fame
)