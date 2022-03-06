cmake_minimum_required (VERSION 2.4)

project (hrotor)

set (CMAKE_CXX_FLAGS "-std=c++14 -Wall -O2 ${CMAKE_CXX_FLAGS}")

aux_source_directory (./ DIR_SRCS)

include_directories (

${PROJECT_SOURCE_DIR}/../include
${PROJECT_SOURCE_DIR})

add_executable (hRotor.x ${DIR_SRCS})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin)

find_package (LAPACK)

if( LAPACK_FOUND)
    set(lapack_libraries ${LAPACK_LIBRARIES})
endif()

target_link_libraries (hRotor.x ${lapack_libraries})

