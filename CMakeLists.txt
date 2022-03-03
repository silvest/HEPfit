cmake_minimum_required(VERSION 3.1)

project(gslpp VERSION 1.0 DESCRIPTION "C++ wrapper for the GSL library")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "-Wall ${GSL_CFLAGS}  -DBOOST_NO_CXX11_RVALUE_REFERENCES")

include_directories(${BOOST_INC})

file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(gslpp ${LIBTYPE} ${srcs})

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${headers} DESTINATION include/RGESolver COMPONENT header)

