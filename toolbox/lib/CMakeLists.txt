cmake_minimum_required( VERSION 3.14 )

project( Splines_Mex )

find_package( Matlab REQUIRED )

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )

## COMPONENTS MX_LIBRARY ENG_LIBRARY MEX_COMPILER )

# set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING INTERNAL FORCE)
set( CMAKE_CXX_STANDARD 11  )

set( SOURCES )
file( GLOB S ${CMAKE_CURRENT_SOURCE_DIR}/../src/*.cc )
foreach(F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../src ${CMAKE_CURRENT_SOURCE_DIR}/../src/Utils )

if ( UNIX OR LINUX )
  if ( NOT APPLE )
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  endif()
else()
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mwindows")
  add_definitions(-mwindows)
endif()

set( MEX_LIST BaseHermite Spline1D Spline2D SplineSet SplineVec )

foreach( BASE ${MEX_LIST})
  matlab_add_mex( NAME ${BASE}MexWrapper SRC ${CMAKE_CURRENT_SOURCE_DIR}/../src_mex/mex_${BASE}MexWrapper.cc ${SOURCES} )

  add_custom_command(
    TARGET ${BASE}MexWrapper
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
    "$<TARGET_FILE:${BASE}MexWrapper>"
    "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
  )

  if ( NOT UNIX )
    target_link_libraries( ${BASE}MexWrapper ws2_32.lib iphlpapi.lib iphlpapi.lib kernel32.lib winpthread )
  endif()

endforeach()
