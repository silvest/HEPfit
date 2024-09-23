cmake_minimum_required(VERSION 3.5)

project(Analysis)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "-Wall ${GSL_CFLAGS} ${ROOT_CFLAGS} -I${BAT_INCLUDE_PATH} ")

foreach(_project ${PROJECTLIST})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../${_project}/src)
endforeach(_project)
include_directories(${BAT_INCLUDE_PATH})
include_directories(${FH_INC})
if(LOOPTOOLS)
  include_directories(${LOOPTOOLS_DIR}/include)
#   set(EXAMPLE_LIBS "-L${LOOPTOOLS_LIB_PATH} -looptools -l${LFORTRAN} ${EXAMPLE_LIBS}")
#   set(EXAMPLE_RPATHS "-Wl,-rpath,${LOOPTOOLS_LIB_PATH} ${EXAMPLE_RPATHS}")
endif(LOOPTOOLS)
include_directories(${BOOST_INC})

set(CMAKE_CXX_LINK_FLAGS "")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${RGESOLVER_LIBS} -L${BAT_LIB_PATH}")

link_directories(${BAT_LIB})
link_directories(${BOOST_LIB_DIR})
if(LOOPTOOLS)
  link_directories(${LOOPTOOLS_LIB_PATH})
endif(LOOPTOOLS)

foreach(_project ${PROJECTLIST})
  set(target_lib_list ${target_lib_list} ${_project})
endforeach(_project)
add_executable(analysis main.cpp)

if(CMAKE_LINK_GROUP_USING_RESCAN_SUPPORTED)
  target_link_libraries(analysis
  "$<LINK_GROUP:RESCAN,${target_lib_list}>" # These libraries will be grouped
  )
else()
  target_link_libraries(analysis ${target_lib_list})
endif()

target_link_libraries(analysis boost_program_options)
target_link_libraries(analysis BATmodels)
target_link_libraries(analysis BAT)
if(LOOPTOOLS)
  target_link_libraries(analysis ooptools)
  target_link_libraries(analysis ${LFORTRAN})
endif(LOOPTOOLS)

if(!${APPLE})
  target_link_libraries(analysis ${LFORTRAN})
endif() 
target_link_libraries(analysis ${ROOT_LDFLAGS})
target_link_libraries(analysis ${ROOT_LIBS})
target_link_libraries(analysis -lMinuit)
target_link_libraries(analysis ${GSL_LIBS})
target_link_libraries(analysis ${RGESOLVER_LIBS})

#add_dependencies(analysis BATBUILD)
#add_dependencies(analysis FHBUILD)
if(MPIBAT)
  target_link_libraries(analysis ${MPI_LIBRARIES})
endif()
if(MPI_COMPILE_FLAGS)
  set_target_properties(analysis PROPERTIES COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(analysis PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
endif()


INSTALL(TARGETS analysis DESTINATION bin COMPONENT executable)

