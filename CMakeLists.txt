cmake_minimum_required(VERSION 2.6)
project(bam_driver)
project(infinite_driver)
project(simulation_driver)
project(simulation_trio)
project(pileup_driver)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Modules")
set(TRIO_MODEL utility.cc read_dependent_data.cc trio_model.cc)
set(EM em_algorithm.cc sufficient_statistics.cc)

include_directories("${CMAKE_SOURCE_DIR}/Eigen/src")
include_directories("${CMAKE_SOURCE_DIR}/bamtools/src")
include_directories("/test")
include_directories("./")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
exec_program(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBRARIES)

add_executable(bam_driver ${TRIO_MODEL} bamtools/src/utils/bamtools_pileup_engine.cpp variant_visitor.cc ${EM} bam_driver.cc)
add_executable(infinite_driver ${TRIO_MODEL} ${EM} infinite_driver.cc)
add_executable(simulation_driver ${TRIO_MODEL} simulation_model.cc simulation_driver.cc)
add_executable(simulation_trio ${TRIO_MODEL} simulation_trio.cc)
add_executable(pileup_driver ${TRIO_MODEL} pileup_utility.cc pileup_driver.cc)

target_link_libraries(bam_driver bamtools bamtools-utils ${GSL_LIBRARIES})
target_link_libraries(infinite_driver ${GSL_LIBRARIES})
target_link_libraries(simulation_driver ${GSL_LIBRARIES})
target_link_libraries(simulation_trio ${GSL_LIBRARIES})
target_link_libraries(pileup_driver ${GSL_LIBRARIES})

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.36.0 COMPONENTS system filesystem unit_test_framework REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

  project(test_utility)
  add_executable(test_utility utility.cc test/test_utility.cc)
  target_link_libraries(test_utility ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_trio_model)
  add_executable(test_trio_model ${TRIO_MODEL} test/test_trio_model.cc)
  target_link_libraries(test_trio_model ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_simulation_model)
  add_executable(test_simulation_model ${TRIO_MODEL} simulation_model.cc test/test_simulation_model.cc)
  target_link_libraries(test_simulation_model ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_pileup_utility)
  add_executable(test_pileup_utility ${TRIO_MODEL} pileup_utility.cc test/test_pileup_utility.cc)
  target_link_libraries(test_pileup_utility ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_em_algorithm)
  add_executable(test_em_algorithm ${TRIO_MODEL} ${EM} test/test_em_algorithm.cc)
  target_link_libraries(test_em_algorithm ${GSL_LIBRARIES} ${Boost_LIBRARIES})

  project(test_sufficient_statistics)
  add_executable(test_sufficient_statistics ${TRIO_MODEL} ${EM} test/test_sufficient_statistics.cc)
  target_link_libraries(test_sufficient_statistics ${GSL_LIBRARIES} ${Boost_LIBRARIES})
endif()