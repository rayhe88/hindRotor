# Set the minimum version of cmake
cmake_minimum_required (VERSION 3.8)

# Set the project name
project (hrotor VERSION "0.0")

execute_process(COMMAND
  git describe --match=NeVeRmAtCh --always --abbrev=7 --dirty
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

configure_file( ${PROJECT_SOURCE_DIR}/../include/version.h.in
                ${PROJECT_SOURCE_DIR}/../include/version.h)
#target_include_directories(hrotor PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                          )


set(CMAKE_CXX_STANDARD 11)

set (CMAKE_CXX_FLAGS "-Werror -Wall -O2 ${CMAKE_CXX_FLAGS}")

aux_source_directory (./ DIR_SRCS)

include_directories (

${PROJECT_SOURCE_DIR}/../include
${PROJECT_SOURCE_DIR})

# Add the executable
add_executable (hRotor.x ${DIR_SRCS})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin)

find_package (LAPACK)

if( LAPACK_FOUND)
    set(lapack_libraries ${LAPACK_LIBRARIES})
endif()

target_link_libraries (hRotor.x ${lapack_libraries})
