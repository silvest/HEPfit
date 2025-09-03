cmake_minimum_required(VERSION 3.10)

project(LoopFunctions)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Wall ${GSL_CFLAGS} ${ROOT_CFLAGS}")
if(LOOPTOOLS)
  set(CMAKE_CXX_FLAGS "-I${LOOPTOOLS_DIR}/include ${CMAKE_CXX_FLAGS}")
endif(LOOPTOOLS)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../gslpp/src)
include_directories(${BOOST_INC})

file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(LoopFunctions ${LIBTYPE} ${srcs})

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${headers} DESTINATION include/HEPfit COMPONENT header)

