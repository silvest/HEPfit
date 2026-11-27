# Writes STAMP only when the content of the installed RGESolver headers changes.
#
# Installing RGESolver keeps the modification time of its source headers, which can be
# older than the HEPfit objects compiled against the previous version, so make and ninja
# would not recompile them. The objects that see the RGESolver headers depend on STAMP.
#
# Usage: cmake -DINC_DIR=<RGESolver include directory> -DSTAMP=<file> -P RGESolverHeaderStamp.cmake

cmake_minimum_required(VERSION 3.10)

# RGESolver.h and the headers it includes with quotes from the same directory, recursively.
set(queue "${INC_DIR}/RGESolver.h")
set(headers "")
while(queue)
  list(GET queue 0 header)
  list(REMOVE_AT queue 0)
  if(NOT EXISTS "${header}" OR header IN_LIST headers)
    continue()
  endif()
  list(APPEND headers "${header}")
  file(STRINGS "${header}" includes REGEX "^[ \t]*#[ \t]*include[ \t]*\"[^\"]+\"")
  foreach(line IN LISTS includes)
    string(REGEX REPLACE "^[ \t]*#[ \t]*include[ \t]*\"([^\"]+)\".*" "\\1" name "${line}")
    list(APPEND queue "${INC_DIR}/${name}")
  endforeach()
endwhile()

set(digest "")
foreach(header IN LISTS headers)
  file(MD5 "${header}" md5)
  string(APPEND digest "${md5}  ${header}\n")
endforeach()
if(NOT headers)
  set(digest "no RGESolver.h in ${INC_DIR}\n")
endif()

set(old "")
if(EXISTS "${STAMP}")
  file(READ "${STAMP}" old)
endif()
if(NOT digest STREQUAL old)
  file(WRITE "${STAMP}" "${digest}")
endif()
