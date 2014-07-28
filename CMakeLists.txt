cmake_minimum_required(VERSION 2.6)
project(bam_driver)
project(simulation_driver)
project(simulation_trio)
project(pileup_driver)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Modules")
find_package(GSL REQUIRED)

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost 1.36.0 COMPONENTS system filesystem unit_test_framework REQUIRED)

if(Boost_FOUND)
  project(test_utility)
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(test_utility utility.cc test_utility.cc)
  target_link_libraries(test_utility ${Boost_LIBRARIES})
endif()

include_directories("${CMAKE_SOURCE_DIR}/Eigen/src")
include_directories("${CMAKE_SOURCE_DIR}/bamtools/src")
include_directories(${GSL_INCLUDE_DIR})
include_directories("./")

EXEC_PROGRAM(${GSL_CONFIG}
          ARGS --libs
          OUTPUT_VARIABLE GSL_LIBRARIES)

add_executable(bam_driver utility.cc read_dependent_data.cc trio_model.cc bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc em_algorithm.cc sufficient_statistics.cc bam_driver.cc)
add_executable(simulation_driver utility.cc read_dependent_data.cc trio_model.cc simulation_model.cc simulation_driver.cc)
add_executable(simulation_trio utility.cc read_dependent_data.cc trio_model.cc simulation_trio.cc)
add_executable(pileup_driver utility.cc read_dependent_data.cc trio_model.cc pileup_utility.cc pileup_driver.cc)

target_link_libraries(bam_driver bamtools bamtools-utils ${GSL_LIBRARIES})
target_link_libraries(simulation_driver ${GSL_LIBRARIES})
target_link_libraries(simulation_trio ${GSL_LIBRARIES})
target_link_libraries(pileup_driver ${GSL_LIBRARIES})