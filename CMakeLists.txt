cmake_minimum_required(VERSION 2.6)
project(bam_driver)

set(CMAKE_CXX_FLAGS "-std=c++11")
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Modules")
find_package( GSL REQUIRED )

include_directories("${CMAKE_SOURCE_DIR}/Eigen/src")
include_directories("${CMAKE_SOURCE_DIR}/bamtools/src")
include_directories(${GSL_INCLUDE_DIR})
include_directories("./")

EXEC_PROGRAM(${GSL_CONFIG}
          ARGS --libs
          OUTPUT_VARIABLE GSL_LIBRARIES)

add_executable(bam_driver utility.cc read_dependent_data.cc trio_model.cc bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc bam_driver.cc)
target_link_libraries(bam_driver bamtools bamtools-utils ${GSL_LIBRARIES})