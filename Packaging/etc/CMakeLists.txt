#--------------------------------------------------------------------
# List of available command-line options:
#   -DDEBUG_MODE:BOOL=<ON or OFF>
#   -DLOCAL_INSTALL_ALL:BOOL=<ON or OFF>
#   -DCMAKE_INSTALL_PREFIX:PATH=<HEPfit installation directory>
#   -DNOMCMC:BOOL=<ON or OFF>
#   -DBOOST_INCLUDE_DIR:PATH=<include path>/boost/
#   -DGSL_CONFIG_DIR:PATH=<gsl-config directory>
#   -DROOT_CONFIG_DIR:PATH=<root-config directory>
#   -DMPIBAT:BOOL=<ON or OFF>
#   -DBAT_INSTALL:BOOL=<ON or OFF>
#   -DBAT_INSTALL_DIR:PATH=<BAT installation directory>
#   -DLOOPTOOLS:BOOL=<ON or OFF>
#   -DLOOPTOOLS_DIR:PATH=<LoopTools directory>
#   -DINTEL_FORTRAN:BOOL=<ON or OFF>
#--------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6)
if(${CMAKE_VERSION} VERSION_LESS 2.8.8)
  set(LIBTYPE STATIC)
else()
  set(LIBTYPE OBJECT)
  #cmake_policy(SET CMP0015 NEW)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(HEPfit)

# default boolean flags
option(DEBUG_MODE "Debug mode" OFF)
option(LOCAL_INSTALL_ALL "Local installation" OFF)
option(NOMCMC "No MCMC mode" OFF)
option(MPIBAT "Use MPI BAT" OFF)
option(BAT_INSTALL "Download and install BAT" OFF)
option(LOOPTOOLS "Use LoopTools library" OFF)
option(INTEL_FORTRAN "Use ifort libraries" OFF)

# Local installation option
if(LOCAL_INSTALL_ALL)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/HEPfit
      CACHE PATH "HEPfit installation directory" FORCE)
  set(BAT_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/BAT
      CACHE PATH "BAT installation directory" FORCE)
  #set(MPIBAT ON CACHE BOOL "Use MPI BAT" FORCE)
  set(BAT_INSTALL ON CACHE BOOL "Download and install BAT" FORCE)
endif()

# ExternalProject module is available since 2.8.0
if(BAT_INSTALL AND ${CMAKE_VERSION} VERSION_LESS 2.8.0)
    message(FATAL_ERROR "The option -DBAT_INSTALL=ON is not applicable for CMake v${CMAKE_VERSION}.")
endif()

# HEPfit installation path
if(NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}
       CACHE PATH "HEPfit installation directory" FORCE)
endif()
string(REGEX REPLACE "^~" "$ENV{HOME}" TMP_DIR ${CMAKE_INSTALL_PREFIX})
set(CMAKE_INSTALL_PREFIX ${TMP_DIR} CACHE PATH "HEPfit installation directory" FORCE)
unset(TMP_DIR)

# BAT installation path
if(BAT_INSTALL_DIR)
  set(BAT_INSTALL_DIR ${BAT_INSTALL_DIR} 
      CACHE PATH "BAT installation directory" FORCE)
else()
  set(BAT_INSTALL_DIR /usr/local 
      CACHE PATH "BAT installation directory")
endif()
if(NOT IS_ABSOLUTE ${BAT_INSTALL_DIR})
   set(BAT_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/${BAT_INSTALL_DIR}
       CACHE PATH "BAT installation directory" FORCE)
endif()
string(REGEX REPLACE "^~" "$ENV{HOME}" TMP_DIR ${BAT_INSTALL_DIR})
set(BAT_INSTALL_DIR ${TMP_DIR} CACHE PATH "BAT installation directory" FORCE)
unset(TMP_DIR)

# Fortran libraries 
if(INTEL_FORTRAN)
  set(LFORTRAN "ifport;ifcore;imf;svml;irc")
  #set(LFORTRAN "ifport;ifcore;imf;svml;irc;ipgo;irc_s")
else()
  set(LFORTRAN "gfortran")
endif()

##########  Colored output  ##########

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColReset    "${Esc}[m")
  set(Red         "${Esc}[31m")
  set(RedBold     "${Esc}[1;31m")
  set(Magenta     "${Esc}[35m")
  set(MagentaBold "${Esc}[1;35m")
endif()

##########  BOOST  ##########

# find_package(Boost) is not used.

set(BOOST_ERROR_MESSAGE "Boost headers could not be found. Please specify the Boost include path with ${RedBold}-DBOOST_INCLUDE_DIR=<include path>/boost/${ColReset}.")

if(NOT BOOST_INCLUDE_DIR)
  find_path(BOOST_INCLUDE_DIR "boost/tokenizer.hpp" DOC "Boost include path")
  if(BOOST_INCLUDE_DIR)
    set(BOOST_INCLUDE_DIR "${BOOST_INCLUDE_DIR}/boost/" CACHE PATH "Boost include path" FORCE)
  else()
    message(FATAL_ERROR "${BOOST_ERROR_MESSAGE}")
  endif()
else()
  set(BOOST_INCLUDE_DIR ${BOOST_INCLUDE_DIR} CACHE PATH "Boost include path" FORCE)
  find_path(TMP_BOOST_INC "tokenizer.hpp" PATHS ${BOOST_INCLUDE_DIR} NO_DEFAULT_PATH)
  if(NOT TMP_BOOST_INC)
    unset(TMP_BOOST_INC CACHE)
    message(FATAL_ERROR "${BOOST_ERROR_MESSAGE}")
  else()
    unset(TMP_BOOST_INC CACHE)
  endif()
endif()

message(STATUS "Found Boost headers in ${BOOST_INCLUDE_DIR}")
string(REGEX REPLACE "/boost$" "" BOOST_INC ${BOOST_INCLUDE_DIR})

##########  GSL  ##########

if(NOT GSL_CONFIG_DIR)
  find_program(GSL_CONFIG gsl-config)
  if(GSL_CONFIG)
    string(REGEX REPLACE "/gsl-config$" "" TMP_PATH ${GSL_CONFIG})
    set(GSL_CONFIG_DIR ${TMP_PATH} CACHE PATH "Path to gsl-config")
    unset(TMP_PATH)
  endif()
else()
  set(GSL_CONFIG_DIR ${GSL_CONFIG_DIR} CACHE PATH "Path to gsl-config" FORCE)
  find_program(GSL_CONFIG gsl-config PATHS ${GSL_CONFIG_DIR} NO_DEFAULT_PATH)
endif()
if(NOT GSL_CONFIG)
  unset(GSL_CONFIG CACHE)
  message(FATAL_ERROR "gsl-confg could not be found. Please specify the directory of gsl-config with ${RedBold}-DGSL_CONFIG_DIR=<directory>${ColReset}.")
else()  
  message(STATUS "Found gsl-config in ${GSL_CONFIG_DIR}")
  execute_process(COMMAND ${GSL_CONFIG} --version
                  OUTPUT_VARIABLE GSL_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GSL version: ${GSL_VERSION}")
  execute_process(COMMAND ${GSL_CONFIG} --cflags
                  OUTPUT_VARIABLE GSL_CFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GSL cflags: ${GSL_CFLAGS}")
  execute_process(COMMAND ${GSL_CONFIG} --libs
                  OUTPUT_VARIABLE GSL_LIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GSL libs: ${GSL_LIBS}")
  string(REGEX MATCH "^-L[^ ]+" GSL_LIB_PATH ${GSL_LIBS})
  if(NOT GSL_LIB_PATH)
    string(REGEX REPLACE "bin$" "lib" GSL_LIB_PATH ${GSL_CONFIG_DIR})
  else()
    string(REGEX REPLACE "^-L" "" GSL_LIB_PATH ${GSL_LIB_PATH})
  endif()
  unset(GSL_CONFIG CACHE)
endif()

##########  ROOT  ##########

if(NOT ROOT_CONFIG_DIR)
  find_program(ROOT_CONFIG root-config DOC "root-config")
  if(ROOT_CONFIG)
    string(REGEX REPLACE "/root-config$" "" TMP_PATH ${ROOT_CONFIG})
    set(ROOT_CONFIG_DIR ${TMP_PATH} CACHE PATH "Path to root-config")
    unset(TMP_PATH)
  endif()
else()
  set(ROOT_CONFIG_DIR ${ROOT_CONFIG_DIR} CACHE PATH "Path to root-config" FORCE)
  find_program(ROOT_CONFIG root-config PATHS ${ROOT_CONFIG_DIR} NO_DEFAULT_PATH DOC "root-config")
endif()
if(NOT ROOT_CONFIG)
  unset(ROOT_CONFIG CACHE)
  message(FATAL_ERROR "root-config could not be found. Please specify the directory of root-config with ${RedBold}-DROOT_CONFIG_DIR=<directory>${ColReset}.")
else()
  message(STATUS "Found root-config in ${ROOT_CONFIG_DIR}")
  execute_process(COMMAND ${ROOT_CONFIG} --version
                  OUTPUT_VARIABLE ROOT_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT version: ${ROOT_VERSION}")
  execute_process(COMMAND ${ROOT_CONFIG} --cflags
                  OUTPUT_VARIABLE ROOT_CFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT cflags: ${ROOT_CFLAGS}")
  execute_process(COMMAND ${ROOT_CONFIG} --libs
                  OUTPUT_VARIABLE ROOT_LIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT libs: ${ROOT_LIBS}")
  execute_process(COMMAND ${ROOT_CONFIG} --ldflags
                  OUTPUT_VARIABLE ROOT_LDFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT ldflags: ${ROOT_LDFLAGS}")
  execute_process(COMMAND ${ROOT_CONFIG} --libdir
                  OUTPUT_VARIABLE ROOT_LIB_PATH
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "/$" "" ROOT_SYS ${ROOT_CONFIG_DIR})
  string(REGEX REPLACE "/bin$" "" ROOT_SYS ${ROOT_SYS})
  unset(ROOT_CONFIG CACHE)
endif()

##########  BAT  ##########

if(NOT NOMCMC)

if(MPIBAT)
  add_definitions(-D_MPI)
  if(CMAKE_CXX_COMPILER MATCHES "(^.*)mpic..$")
    set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "MPI CXX compiler wrapper" FORCE)
  endif()
  if(${CMAKE_VERSION} VERSION_LESS 2.8.5)
    # MPI_CXX_COMPILER is available since 2.8.5
    set(MPI_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "MPI CXX compiler wrapper" FORCE)
  endif()
  # FindMPI attempts to determine MPI variables using MPI_CXX_COMPILER (or MPI_COMPILER) if available. 
  find_package(MPI)
  if(${CMAKE_VERSION} VERSION_LESS 2.8.5)
      # MPI_CXX_FOUND, etc. are available since 2.8.5
      set(MPI_CXX_FOUND ${MPI_FOUND})
      set(MPI_CXX_COMPILER ${MPI_COMPILER} CACHE FILEPATH "MPI CXX compiler wrapper" FORCE)
      set(MPI_CXX_INCLUDE_PATH ${MPI_INCLUDE_PATH} CACHE PATH "MPI CXX include path" FORCE)
  endif()
  if(NOT MPI_CXX_FOUND)
    string(REGEX REPLACE "/bin/mpic..$" "" TMP_MPI_PATH ${MPI_CXX_COMPILER}) 
    find_path(TMP_MPI_BIN_PATH "mpicxx" PATHS ${TMP_MPI_PATH}/bin)
    find_path(TMP_MPI_INC_PATH "mpi.h" PATHS ${TMP_MPI_PATH}/include)
    unset(TMP_MPI_PATH)
    if(TMP_MPI_BIN_PATH AND TMP_MPI_INC_PATH)
      message(STATUS "Found MPI binaries in ${TMP_MPI_BIN_PATH} and headers in ${TMP_MPI_INC_PATH}")
      set(MPI_CXX_COMPILER ${TMP_MPI_BIN_PATH}/mpicxx CACHE FILEPATH "MPI CXX compiler wrapper" FORCE)
      set(MPI_CXX_INCLUDE_PATH ${TMP_MPI_INC_PATH} CACHE PATH "MPI CXX include path" FORCE)
      include_directories(${MPI_CXX_INCLUDE_PATH})
      unset(TMP_MPI_BIN_PATH CACHE)
      unset(TMP_MPI_INC_PATH CACHE)
    else()
      unset(TMP_MPI_BIN_PATH CACHE)
      unset(TMP_MPI_INC_PATH CACHE)
      message(FATAL_ERROR "MPI support has been requested but implementation could not be found. Please specify the location of the mpicxx compiler with ${RedBold}-DMPI_CXX_COMPILER=<path to mpi>/mpicxx${ColReset} or ${RedBold}-DCMAKE_CXX_COMPILER=<path to mpi>/mpicxx${ColReset}.")
    endif()
  endif(NOT MPI_CXX_FOUND)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  string(REGEX REPLACE "/bin/mpic..$" "" TMP_MPI_PATH ${MPI_CXX_COMPILER}) 
  if(IS_DIRECTORY "${TMP_MPI_PATH}/lib")
    set(MPI_LIB_PATH ${TMP_MPI_PATH}/lib)
  elseif(IS_DIRECTORY "${TMP_MPI_PATH}/lib64")
    set(MPI_LIB_PATH ${TMP_MPI_PATH}/lib64)
  else()
    message(STATUS "Not defined MPI_LIB_PATH")
  endif()
  unset(TMP_MPI_PATH)
endif(MPIBAT)

message(STATUS "Looking for BAT in ${BAT_INSTALL_DIR}")
find_path(IBAT "BCModel.h" PATHS ${BAT_INSTALL_DIR}/include/BAT NO_DEFAULT_PATH)
find_library(LBAT BAT PATHS ${BAT_INSTALL_DIR} NO_DEFAULT_PATH PATH_SUFFIXES lib lib32 lib64)
if(IBAT AND LBAT)
  string(REGEX REPLACE "/BAT$" "" BAT_INCLUDE_PATH ${IBAT})
  string(REGEX REPLACE "/libBAT.*$" "" BAT_LIB_PATH ${LBAT})
  message(STATUS "Found BAT headers in ${BAT_INCLUDE_PATH} and libs in ${BAT_LIB_PATH}")
  if(BAT_INSTALL)
    message(STATUS "${RedBold}BAT will not be installed, since it's already installed in ${BAT_INSTALL_DIR}.${ColReset}")
  endif()
  if(MPIBAT)
    message(STATUS "${RedBold}Please make sure that the installed BAT library has been compiled with our patch for MPI support.${ColReset}")
  endif()
else(IBAT AND LBAT)
  if(NOT BAT_INSTALL)
    unset(IBAT CACHE)
    unset(LBAT CACHE)
    message(FATAL_ERROR "BAT could not be found in ${BAT_INSTALL_DIR}. Please supply the correct path to BAT with the option ${RedBold}-DBAT_INSTALL_DIR=<BAT installation path>${ColReset}. Alternatively, setting the option ${RedBold}-DBAT_INSTALL=ON${ColReset}, I'll download and install BAT at the default location /usr/local/, where you can specify a custom installatation location with ${RedBold}-DBAT_INSTALL_DIR=<BAT custom installation path>${ColReset}.")  
  endif()
endif(IBAT AND LBAT)

if(BAT_INSTALL)
  if(NOT IBAT OR NOT LBAT)
    message(STATUS "Not found BAT in ${BAT_INSTALL_DIR}")
    message(STATUS "${RedBold}I'll download and install BAT into ${BAT_INSTALL_DIR}.${ColReset}")
    include(ExternalProject)
    set(BAT_VERSION "0.9.4.1")
    set(BAT_FILENAME "BAT-${BAT_VERSION}.tar.gz")
    find_file(BAT_URL ${BAT_FILENAME} PATHS ${CMAKE_CURRENT_BINARY_DIR} NO_DEFAULT_PATH)
    if(NOT BAT_URL)
      set(BAT_URL https://github.com/bat/bat/releases/download/v${BAT_VERSION}/${BAT_FILENAME})
    else()
      message(STATUS "Found ${BAT_FILENAME} in the current directory")
    endif()
    set(BAT_MD5 e1a09dc5e379dd5b5d665729ea45b535)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/BAT_make_wrapper.sh.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/BAT_make_wrapper.sh" @ONLY)
    # The make command to the top-level make will be inherited to the BAT make. 
    set(BAT_MAKE $(MAKE))
    set(BAT_CONFIGURE "./configure")
    if(BAT_INSTALL_DIR)
      set(BAT_CONFIGURE "${BAT_CONFIGURE} --prefix=${BAT_INSTALL_DIR}")
    endif()
    if(ROOT_SYS)
      set(BAT_CONFIGURE "${BAT_CONFIGURE} --with-rootsys=${ROOT_SYS}")
    endif()
    if(MPIBAT)
      message(STATUS "${RedBold}Enable MPI support${ColReset}")
      set(PATCH patch -p1 -N < ${CMAKE_CURRENT_SOURCE_DIR}/BAT_mpi_patch.txt)
      ExternalProject_Add(
        BATBUILD
        PREFIX BATBUILD
        TMP_DIR BATBUILD/tmp
        STAMP_DIR BATBUILD/stamp
        SOURCE_DIR BATBUILD/src
        BINARY_DIR BATBUILD/src
        INSTALL_DIR BATBUILD/src
        DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
        URL ${BAT_URL}
        URL_MD5 ${BAT_MD5}
        PATCH_COMMAND ${PATCH}
        CONFIGURE_COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/BAT_make_wrapper.sh ${BAT_CONFIGURE} ${MPI_CXX_COMPILER}
        BUILD_COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/BAT_make_wrapper.sh ${BAT_MAKE}
        INSTALL_COMMAND make install
        )
    else(MPIBAT)
      ExternalProject_Add(
        BATBUILD
        PREFIX BATBUILD
        TMP_DIR BATBUILD/tmp
        STAMP_DIR BATBUILD/stamp
        SOURCE_DIR BATBUILD/src
        BINARY_DIR BATBUILD/src
        INSTALL_DIR BATBUILD/src
        DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
        URL ${BAT_URL}
        URL_MD5 ${BAT_MD5}
        CONFIGURE_COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/BAT_make_wrapper.sh ${BAT_CONFIGURE}
        BUILD_COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/BAT_make_wrapper.sh ${BAT_MAKE}
        INSTALL_COMMAND make install
        )
    endif(MPIBAT)
    set(BAT_INCLUDE_PATH ${BAT_INSTALL_DIR}/include)
    set(BAT_LIB_PATH ${BAT_INSTALL_DIR}/lib)
  endif(NOT IBAT OR NOT LBAT)
endif(BAT_INSTALL)

unset(IBAT CACHE)
unset(LBAT CACHE)

endif(NOT NOMCMC)

##########  LoopTools  ##########

if(LOOPTOOLS)
  add_definitions(-DUSE_LOOPTOOLS)
  if(NOT LOOPTOOLS_DIR)
    find_path(TMP_INC "clooptools.h")
    find_library(TMP_LIB ooptools)
  else()
    set(LOOPTOOLS_DIR ${LOOPTOOLS_DIR} CACHE PATH "LoopTools directory" FORCE)
    find_path(TMP_INC "clooptools.h" PATHS ${LOOPTOOLS_DIR}/include NO_DEFAULT_PATH)
    find_library(TMP_LIB ooptools PATHS ${LOOPTOOLS_DIR} NO_DEFAULT_PATH PATH_SUFFIXES lib lib32 lib64)
  endif()
  if(TMP_INC AND TMP_LIB)
    if(NOT LOOPTOOLS_DIR)
      string(REGEX REPLACE "/include$" "" TMP_INC_PATH ${TMP_INC})
      set(LOOPTOOLS_DIR "${TMP_INC_PATH}" CACHE PATH "LoopTools directory")
    endif()
    string(REGEX REPLACE "/libooptools.*$" "" TMP_LIB_PATH ${TMP_LIB})
    set(LOOPTOOLS_LIB_PATH ${TMP_LIB_PATH})
    unset(TMP_LIB_PATH)
    message(STATUS "Found LoopTools headers in ${LOOPTOOLS_DIR}/include and libs in ${LOOPTOOLS_LIB_PATH}")
  else()
    message(FATAL_ERROR "LoopTools headers and/or libs could not be found. Please specify the LoopTools path with ${RedBold}-DLOOPTOOLS_DIR=<LoopTools installation path>${ColReset}.")
  endif()
  unset(TMP_INC CACHE)
  unset(TMP_LIB CACHE)
endif(LOOPTOOLS)

##########  FeynHiggs  ##########
#
#option(FH_INSTALL "Install FeynHiggs" OFF)
#if(NOT FH_INSTALL_DIR)
#  set(FH_INSTALL_DIR /usr/local)
#endif()
#if(NOT IS_ABSOLUTE ${FH_INSTALL_DIR})
#   set(FH_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/${FH_INSTALL_DIR})
#endif()
#STRING(REGEX REPLACE "^~" "$ENV{HOME}" FH_INSTALL_DIR ${FH_INSTALL_DIR})
#
#message(STATUS "looking for FeynHiggs in ${FH_INSTALL_DIR}")
#find_path(IFH "CFeynHiggs.h" PATHS ${FH_INSTALL_DIR}/include NO_DEFAULT_PATH)
#find_library(LFH NAMES libFH.a FH PATHS ${FH_INSTALL_DIR} NO_DEFAULT_PATH PATH_SUFFIXES lib lib32 lib64)
#if(IFH AND LFH)
#  set(FH_INC ${IFH})
#  string(REGEX REPLACE "/libFH.*$" "" FH_LIB ${LFH})
#  message(STATUS "Found FeynHiggs include in ${FH_INC}")
#  message("   and libs in ${FH_LIB}")
#else(IFH AND LFH)
#  if(NOT FH_INSTALL)
#    message(FATAL_ERROR "FeynHiggs include files and libs could not be found. Please specify the correct path to FeynHiggs with -DFH_INSTALL_DIR=<FeynHiggs installation path> or install FeynHiggs with the option -DFH_INSTALL=on, where I'll download and install FeynHiggs at the default location /usr/local/. If you want me to install FeynHiggs at a custom location, you can specify it with -DFH_INSTALL_DIR=<FeynHiggs installation path>")
#  endif()
#endif(IFH AND LFH)
#
#if(FH_INSTALL)
#  if(NOT IFH OR NOT LFH)
#    message(STATUS "${MagentaBold}I'll download and install FeynHiggs into ${FH_INSTALL_DIR}.${ColReset}")
#    set(CONFIGURE_FH ./configure --prefix=${FH_INSTALL_DIR})
#    set(FH_INC ${FH_INSTALL_DIR}/include)
#    set(FH_LIB ${FH_INSTALL_DIR}/lib ${FH_INSTALL_DIR}/lib32 ${FH_INSTALL_DIR}/lib64)
#    include(ExternalProject)
#    ExternalProject_Add(
#      FHBUILD
#      PREFIX FHBUILD
#      SOURCE_DIR FHBUILD/src
#      STAMP_DIR FHBUILD/src
#      BINARY_DIR FHBUILD/src
#      INSTALL_DIR FHBUILD/src
#      URL http://wwwth.mpp.mpg.de/members/heinemey/feynhiggs/newversion/FeynHiggs-2.10.0.tar.gz
#      CONFIGURE_COMMAND ${CONFIGURE_FH}
#      BUILD_COMMAND make
#      INSTALL_COMMAND make install
#      )
#  endif()
#endif()
#
##message("FH_INC: ${FH_INC}")
##message("FH_LIB: ${FH_LIB}")
##message("LFH: ${LFH}")

##########  DEBUG  ##########

if(DEBUG_MODE)
  message("\nCache variables:")
  message("  CMAKE_CXX_COMPILER:FILEPATH = ${CMAKE_CXX_COMPILER}")
  message("  CMAKE_INSTALL_PREFIX:PATH = ${CMAKE_INSTALL_PREFIX}")
  message("  DEBUG_MODE:BOOL = ${DEBUG_MODE}")
  message("  LOCAL_INSTALL_ALL:BOOL = ${LOCAL_INSTALL_ALL}")
  message("  NOMCMC:BOOL = ${NOMCMC}")
  message("  BOOST_INCLUDE_DIR:PATH = ${BOOST_INCLUDE_DIR}")
  message("  GSL_CONFIG_DIR:PATH = ${GSL_CONFIG_DIR}")
  message("  ROOT_CONFIG_DIR:PATH = ${ROOT_CONFIG_DIR}")
  message("  MPI_CXX_COMPILER:FILEPATH = ${MPI_CXX_COMPILER}")
  #message("  MPI_CXX_INCLUDE_PATH:PATH = ${MPI_CXX_INCLUDE_PATH}")
  message("  MPIBAT:BOOL = ${MPIBAT}")
  message("  BAT_INSTALL:BOOL = ${BAT_INSTALL}")
  message("  BAT_INSTALL_DIR:PATH = ${BAT_INSTALL_DIR}")
  message("  LOOPTOOLS:BOOL = ${LOOPTOOLS}")
  message("  LOOPTOOLS_DIR:PATH = ${LOOPTOOLS_DIR}")
  message("  INTEL_FORTRAN:BOOL = ${INTEL_FORTRAN}")

  message("\nNormal variables:")
  message("  LIBTYPE = ${LIBTYPE}")
  message("  BOOST_INC = ${BOOST_INC}")
  message("  GSL_VERSION = ${GSL_VERSION}")
  message("  GSL_CFLAGS = ${GSL_CFLAGS}")
  message("  GSL_LIBS = ${GSL_LIBS}")
  message("  GSL_LIB_PATH = ${GSL_LIB_PATH}")
  message("  ROOT_VERSION = ${ROOT_VERSION}")
  message("  ROOT_CFLAGS = ${ROOT_CFLAGS}")
  message("  ROOT_LIBS = ${ROOT_LIBS}")
  message("  ROOT_LDFLAGS = ${ROOT_LDFLAGS}")
  message("  ROOT_LIB_PATH = ${ROOT_LIB_PATH}")
  message("  ROOT_SYS = ${ROOT_SYS}")
  message("  MPI_LIB_PATH=${MPI_LIB_PATH}")
  message("  BAT_FILENAME = ${BAT_FILENAME}")
  message("  BAT_URL = ${BAT_URL}")
  message("  BAT_MD5 = ${BAT_MD5}")
  message("  BAT_CONFIGURE = ${BAT_CONFIGURE}")
  message("  BAT_MAKE = ${BAT_MAKE}")
  message("  BAT_INCLUDE_PATH = ${BAT_INCLUDE_PATH}")
  message("  BAT_LIB_PATH = ${BAT_LIB_PATH}")
  message("  LOOPTOOLS_LIB_PATH = ${LOOPTOOLS_LIB_PATH}")
  message("  LFORTRAN = ${LFORTRAN}")
  message("")

endif(DEBUG_MODE)

##############################

#set(PROJECTLIST ComputeObservables EventGeneration InputParser EW Flavour LeptonFlavour NewPhysics Observables GeneralSUSY pMSSM MFV SUSY THDM HiggsExtensions StandardModel LoopFunctions SUSYMassInsertion gslpp)
set(PROJECTLIST ComputeObservables InputParser EW Flavour NewPhysics Observables HiggsExtensions StandardModel LoopFunctions gslpp)
if(NOT NOMCMC)
  set(PROJECTLIST MonteCarlo ${PROJECTLIST})
endif(NOT NOMCMC)

# sub-directories
foreach(_project ${PROJECTLIST})
  add_subdirectory(${_project})
endforeach(_project)

# HEPfit library libHEPfit.a
if(${LIBTYPE} STREQUAL "OBJECT") 
  # for CMake v2.8.8 or later
  foreach(_project ${PROJECTLIST})
    set(OBJECTLIST ${OBJECTLIST} $<TARGET_OBJECTS:${_project}>)
  endforeach(_project)
  add_library(HEPfit STATIC ${OBJECTLIST})
else()
  # for older CMake
  set(CMAKE_CXX_FLAGS "-Wall ${ROOT_CFLAGS} ${GSL_CFLAGS}")
  if(NOT NOMCMC)
    include_directories(${BAT_INCLUDE_PATH})
  endif(NOT NOMCMC)
  if(LOOPTOOLS)
    include_directories(${LOOPTOOLS_DIR}/include)
  endif(LOOPTOOLS)
  include_directories(${BOOST_INC})
  foreach(_project ${PROJECTLIST})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${_project}/src)
    file(GLOB SOURCES ${_project}/src/*.cpp)
    set(SRCLIST ${SRCLIST} ${SOURCES})
  endforeach(_project)
  add_library(HEPfit STATIC ${SRCLIST})
endif()
if(BAT_INSTALL AND NOT NOMCMC)
  add_dependencies(HEPfit BATBUILD)
endif(BAT_INSTALL AND NOT NOMCMC)
INSTALL(TARGETS HEPfit DESTINATION lib COMPONENT library)

# HEPfit header file HEPfit.h 
list(REVERSE PROJECTLIST)
foreach(_project ${PROJECTLIST})
  file(GLOB HEADERS ${_project}/src/*.h)
  foreach(_header ${HEADERS})
    get_filename_component(_tmp ${_header} NAME)
    if(NOT ${_tmp} MATCHES "InputParameters.h")
      set(INCLUDE_LIST "${INCLUDE_LIST}\#include\ \"${_tmp}\"\n")
      if(NOT ${_tmp} MATCHES "MonteCarlo.h" AND NOT ${_tmp} MATCHES "MonteCarloEngine.h")
        set(NOMCMC_INCLUDE_LIST "${NOMCMC_INCLUDE_LIST}\#include\ \"${_tmp}\"\n")
      endif()
    endif()
  endforeach(_header)
endforeach(_project)
if(NOT NOMCMC)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/HEPfit.h.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/HEPfit.h" @ONLY)
else(NOT NOMCMC)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/HEPfit_noMCMC.h.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/HEPfit.h" @ONLY)
endif(NOT NOMCMC)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/HEPfit.h 
        DESTINATION include/HEPfit COMPONENT header)
      
# pkg-config (FindPkgConfig in CMake 2.8.8 contains a bug)
#find_package(PkgConfig)
#if(PKG_CONFIG_FOUND)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/HEPfit.pc.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/HEPfit.pc" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/hepfit-config.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/hepfit-config" @ONLY)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/hepfit-config ${CMAKE_CURRENT_BINARY_DIR}/HEPfit.pc 
          DESTINATION bin COMPONENT hepfit-config)
#endif()

# example programs
set(EXAMPLE_CXXFLAGS "-I${BOOST_INC}")
set(EXAMPLE_LIBS " ")
set(EXAMPLE_RPATHS "-Wl,-rpath,${GSL_LIB_PATH}")
set(EXAMPLE_RPATHS "-Wl,-rpath,${ROOT_LIB_PATH} ${EXAMPLE_RPATHS}")
if(MPI_LIB_PATH)
  set(EXAMPLE_RPATHS "-Wl,-rpath,${MPI_LIB_PATH} ${EXAMPLE_RPATHS}")
endif(MPI_LIB_PATH)
if(LOOPTOOLS)
  set(EXAMPLE_CXXFLAGS "-I${LOOPTOOLS_DIR}/include ${EXAMPLE_CXXFLAGS}")
  set(EXAMPLE_LIBS "-L${LOOPTOOLS_LIB_PATH} -looptools -l${LFORTRAN} ${EXAMPLE_LIBS}")
  set(EXAMPLE_RPATHS "-Wl,-rpath,${LOOPTOOLS_LIB_PATH} ${EXAMPLE_RPATHS}")
endif(LOOPTOOLS)
#
#set(EXAMPLE_CXXFLAGS "-I${FH_INC} ${EXAMPLE_CXXFLAGS}")
#set(EXAMPLE_LIBS "-L${FH_LIB} -lFH ${EXAMPLE_LIBS}")
#set(EXAMPLE_RPATHS "-Wl,-rpath,${FH_LIB} ${EXAMPLE_RPATHS}")
#
set(EXAMPLE_NOMCMC_CXXFLAGS "${EXAMPLE_CXXFLAGS}")
set(EXAMPLE_NOMCMC_LIBS "${EXAMPLE_LIBS}")
set(EXAMPLE_NOMCMC_RPATHS "${EXAMPLE_RPATHS}")
if(NOT NOMCMC)
  if(MPIBAT)
    set(EXAMPLE_MCMC_CXX "${MPI_CXX_COMPILER}")
    set(EXAMPLE_CXXFLAGS "-D_MPI ${EXAMPLE_CXXFLAGS}")
  else(MPIBAT)
    set(EXAMPLE_MCMC_CXX "${CMAKE_CXX_COMPILER}")
  endif(MPIBAT)
  set(EXAMPLE_CXXFLAGS "-I${BAT_INCLUDE_PATH} ${EXAMPLE_CXXFLAGS}")
  set(EXAMPLE_LIBS "-L${BAT_LIB_PATH} -lBATmodels -lBAT ${EXAMPLE_LIBS}")
  set(EXAMPLE_RPATHS "-Wl,-rpath,${BAT_LIB_PATH} ${EXAMPLE_RPATHS}")
endif(NOT NOMCMC)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/LibMode_config/Makefile.in"
               "${CMAKE_CURRENT_BINARY_DIR}/examples/LibMode_config/Makefile" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/LibMode_config/libmode_config.cpp"
               "${CMAKE_CURRENT_BINARY_DIR}/examples/LibMode_config/libmode_config.cpp" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/LibMode_header/Makefile.in"
               "${CMAKE_CURRENT_BINARY_DIR}/examples/LibMode_header/Makefile" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/LibMode_header/libmode_header.cpp"
               "${CMAKE_CURRENT_BINARY_DIR}/examples/LibMode_header/libmode_header.cpp" @ONLY)
file(GLOB CONF_FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples-src/config/*.conf)
foreach(_confFile ${CONF_FILES})
  get_filename_component(_confFileName ${_confFile} NAME)
  configure_file("${_confFile}" 
                 "${CMAKE_CURRENT_BINARY_DIR}/examples/config/${_confFileName}" COPYONLY)
endforeach(_confFile)
file(GLOB INPUT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples-src/config/input/*.root)
foreach(_inputFile ${INPUT_FILES})
  get_filename_component(_inputFileName ${_inputFile} NAME)
  configure_file("${_inputFile}" 
                 "${CMAKE_CURRENT_BINARY_DIR}/examples/config/input/${_inputFileName}" COPYONLY)
endforeach(_inputFile)

if(NOT NOMCMC)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/MonteCarloMode/Makefile.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/examples/MonteCarloMode/Makefile" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/MonteCarloMode/MCMC.cpp"
                 "${CMAKE_CURRENT_BINARY_DIR}/examples/MonteCarloMode/MCMC.cpp" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/MonteCarloMode/MonteCarlo.conf"
                 "${CMAKE_CURRENT_BINARY_DIR}/examples/MonteCarloMode/MonteCarlo.conf" COPYONLY)

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/examples-src/myModel/Makefile.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/examples/myModel/Makefile" @ONLY)
  set(_myDirs . config src)
  foreach(_myDir ${_myDirs})
    set(_myExts h cpp conf)
    foreach(_myExt ${_myExts})
      file(GLOB _myFiles ${CMAKE_CURRENT_SOURCE_DIR}/examples-src/myModel/${_myDir}/*.${_myExt})
      foreach(_myFile ${_myFiles})
        get_filename_component(_myFileName ${_myFile} NAME)
        configure_file("${_myFile}" 
                       "${CMAKE_CURRENT_BINARY_DIR}/examples/myModel/${_myDir}/${_myFileName}" COPYONLY)
        #message(STATUS "${_myFile}")
      endforeach(_myFile)
    endforeach(_myExt)
  endforeach(_myDir)
endif(NOT NOMCMC)
message(STATUS "Generated example codes in examples/")

message(STATUS "CMAKE_INSTALL_PREFIX has been set to ${CMAKE_INSTALL_PREFIX}")

# make uninstall
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} 
                  -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

