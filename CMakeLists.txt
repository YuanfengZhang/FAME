cmake_minimum_required(VERSION 3.12)
project(FAME VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags for optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find required dependencies
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

# Optional: zstd support for .zst compressed files
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
	pkg_check_modules(ZSTD libzstd)
endif()

# Fallback: try to find zstd manually
if(NOT ZSTD_FOUND)
	find_path(ZSTD_INCLUDE_DIR zstd.h)
	find_library(ZSTD_LIBRARY NAMES zstd libzstd)

	if(ZSTD_INCLUDE_DIR AND ZSTD_LIBRARY)
		set(ZSTD_FOUND TRUE)
		set(ZSTD_LIBRARIES ${ZSTD_LIBRARY})
		set(ZSTD_INCLUDE_DIRS ${ZSTD_INCLUDE_DIR})
	endif()
endif()

# Report zstd status
if(ZSTD_FOUND)
	message(STATUS "zstd found: ${ZSTD_LIBRARIES}")
	message(STATUS "zstd include: ${ZSTD_INCLUDE_DIRS}")
	set(ZSTD_SUPPORT ON)
else()
	message(WARNING "zstd not found, .zst support disabled. Install libzstd-dev for zstd support.")
	set(ZSTD_SUPPORT OFF)
endif()

# Source files
set(SOURCES
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

# Add zstdstream if zstd is available
if(ZSTD_SUPPORT)
	list(APPEND SOURCES zstdstream/zstdstream.C)
endif()

# Main executable
add_executable(fame ${SOURCES})

target_include_directories(fame PRIVATE
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/sparsehash/include/usr/local/include
	${CMAKE_SOURCE_DIR}/hopscotch-map/include/tsl
	${CMAKE_SOURCE_DIR}/spaced_nthash
)

# Link required libraries
target_link_libraries(fame PRIVATE
	OpenMP::OpenMP_CXX
	ZLIB::ZLIB
)

# Link zstd if available
if(ZSTD_SUPPORT)
	target_include_directories(fame PRIVATE ${ZSTD_INCLUDE_DIRS})
	target_link_libraries(fame PRIVATE ${ZSTD_LIBRARIES})
	target_compile_definitions(fame PRIVATE ZSTD_SUPPORT)
endif()

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