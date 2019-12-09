############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.5)
SET( CMAKE_VERBOSE_MAKEFILE OFF )
SET( CMAKE_INSTALL_MESSAGE NEVER )

IF ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  SET( DEBUG_SUFFIX "_debug" )
ELSE()
  SET( DEBUG_SUFFIX )
ENDIF()

IF ( BUILD_SHARED )
  IF ( APPLE )
    SET( TARGET Splines_osx${DEBUG_SUFFIX} )
  ELSEIF( LINUX OR UNIX )
    SET( TARGET Splines_linux${DEBUG_SUFFIX} )
  ELSE()
    SET( TARGET Splines_win_${BITS}${DEBUG_SUFFIX} )
    SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
  ENDIF()
ENDIF()

IF ( APPLE )
  SET( TARGETS Splines_osx_static${DEBUG_SUFFIX} )
ELSEIF( LINUX OR UNIX )
  SET( TARGETS Splines_linux_static${DEBUG_SUFFIX} )
ELSE()
  SET( TARGETS Splines_win_${BITS}_static${DEBUG_SUFFIX} )
ENDIF()

PROJECT( Splines CXX C )

INCLUDE( ./CMakeLists-cflags.txt )
MESSAGE( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")

SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})
SET( HEADERS )
FILE( GLOB S ./src/*.hh )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})
FILE( GLOB S ./src/*.h )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

IF( NOT DEFINED GC_DIR )
  SET( GC_DIR "./GC" )
ENDIF()

FILE( GLOB GC_DIR_GLOB ${GC_DIR} )
MESSAGE( STATUS "GC_DIR_GLOB = ${GC_DIR_GLOB}" )

INCLUDE_DIRECTORIES( src ${GC_DIR_GLOB}/lib/include )
LINK_DIRECTORIES( ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll )
IF ( APPLE )
  FIND_LIBRARY(
    GC_LIB NAMES GenericContainer_osx_static${DEBUG_SUFFIX}
    HINTS ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll
  )
ELSEIF( LINUX OR UNIX )
  FIND_LIBRARY(
    GC_LIB NAMES GenericContainer_linux_static${DEBUG_SUFFIX}
    HINTS ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll
  )
ELSE()
  FIND_LIBRARY(
    GC_LIB NAMES GenericContainer_win_${BITS}_static${DEBUG_SUFFIX}
    HINTS ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll
  )
ENDIF()

IF ( BUILD_SHARED )
  IF ( APPLE )
    FIND_LIBRARY(
      GC_LIB_SHARED NAMES GenericContainer_osx${DEBUG_SUFFIX}
      HINTS ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll
    )
  ELSEIF( LINUX OR UNIX )
    FIND_LIBRARY(
      GC_LIB_SHARED NAMES GenericContainer_linux${DEBUG_SUFFIX}
      HINTS ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll
    )
  ELSE()
    FIND_LIBRARY(
      GC_LIB_SHARED NAMES GenericContainer_win_${BITS}${DEBUG_SUFFIX}
      HINTS ${GC_DIR_GLOB}/lib/lib ${GC_DIR_GLOB}/lib/dll
    )
  ENDIF()
ENDIF()

SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib )

IF ( BUILD_SHARED )
  ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
  ADD_LIBRARY( ${TARGET}  SHARED ${SOURCES} ${HEADERS} )
  TARGET_LINK_LIBRARIES( ${TARGET} ${GC_LIB_SHARED} )
  INSTALL(
    TARGETS ${TARGET} ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ELSE()
  ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
  INSTALL(
    TARGETS ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ENDIF()

INSTALL( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include )

IF( BUILD_EXECUTABLE )

  SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FILE( MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )

  FILE( GLOB S ./tests/*.cc )
  FOREACH (F ${S})
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    get_filename_component( TEXE ${RF} NAME_WE )
    ADD_EXECUTABLE( ${TEXE} ${RF} )
    TARGET_LINK_LIBRARIES( ${TEXE} ${TARGETS} ${GC_LIB} )
  ENDFOREACH (F ${S})
ENDIF()

MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
MESSAGE( STATUS "Sources      = ${SOURCES}" )
MESSAGE( STATUS "Headers      = ${HEADERS}" )
MESSAGE( STATUS "GC_LIB       = ${GC_LIB}" )
MESSAGE( STATUS "GC_DIR       = ${GC_DIR}" )
MESSAGE( STATUS "GC_DIR_GLOB  = ${GC_DIR_GLOB}" )
