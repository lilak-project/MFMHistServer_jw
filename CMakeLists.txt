cmake_minimum_required(VERSION 2.8)
project(MFMHistServer)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/local/etc/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
add_definitions( -DTIXML_USE_STL )

include_directories(${MFMHistServer_SOURCE_DIR}/source)

set(GETDIR "/usr/local/get")
set(GRUDIR "/usr/local/gru")
set(GETDIR2 "/opt/software/get_src/20210208/MultiFrame/src")
set(LIBGRUDIR "${GRUDIR}")

find_package(ROOT REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} -std=c++14 -pthread -Wno-deprecated-declarations")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lncurses" )
include_directories(${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR})
link_directories(${LIBGRUDIR}/lib)
link_directories(${ROOT_LIBRARY_DIR})

# Find FFTW3
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
link_directories(${FFTW_LIBRARIES})

include_directories(${GETDIR}/include)
include_directories(${GRUDIR}/include)
include_directories(${GRUDIR}/vigru)
include_directories(${GRUDIR}/gsoap)

include_directories(${GETDIR2})

link_directories(${GRUDIR}/lib)
link_directories(${GETDIR}/lib)
#include_directories(${GETDIR2})

set(SRCFILES source/HistServer.cpp source/jsoncpp.cpp source/WaveletNew.cpp)
set(HDRFILES source/HistServer.h)

add_executable(run_mfm_converter main.cpp ${HDRFILES} ${SRCFILES})
target_link_libraries(run_mfm_converter ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} Html Spectrum XMLIO GRU MultiFrame)
