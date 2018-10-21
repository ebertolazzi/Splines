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

SET( SPLINES Splines )

PROJECT( ${SPLINES} CXX C )

FILE( GLOB S ./CMakeLists-cflags.txt )
IF(EXISTS ${S})
  INCLUDE( ${S} )
ELSE()
  SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
  SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
  SET( CMAKE_BUILD_TYPE                  Release )
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
    SET( CMAKE_CXX_FLAGS "-std=c++11 " )
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Wall -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -funroll-loops -O3 -g0 " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -Wall -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -O0 -g3 -DMECHATRONIX_DEBUG " )
  ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
    SET( CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ " )
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -funroll-loops -O3 -g0 -fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation -O0 -gfull -DMECHATRONIX_DEBUG -fPIC -Weverything -Wno-weak-vtables -Wno-implicit-fallthrough -Wno-documentation-unknown-command -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-reserved-id-macro -Wno-documentation " )
  ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    SET( CMAKE_CXX_FLAGS_RELEASE "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC /O2 /MD " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS /D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC /Od /Ob0 /MDd /Zi /RTC1 /D_DEBUG " )
  ELSE()
    MESSAGE( FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
  ENDIF()
ENDIF()

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

FILE( GLOB S ./GC/lib/include )
IF(EXISTS ${S})
  INCLUDE_DIRECTORIES( src ${S} )
ENDIF()

FILE( GLOB S ./GC/lib )
IF(EXISTS ${S})
  LINK_DIRECTORIES( ${S} )
ENDIF()

FIND_LIBRARY( GC_LIB NAMES GenericContainer libGenericContainer HINTS ./GC/lib )

IF (WIN32)
  IF( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    IF ( NOT DEFINED BITS)
      SET( BITS "x64" )
    ENDIF()
    INCLUDE_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include" )
    IF ( ${BITS} STREQUAL "x64" )
      SET( WINXXBITS Win64 )
      LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64" )
    ELSE()
      SET( WINXXBITS Win32 )
      LINK_DIRECTORIES( "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/ia32" )
    ENDIF()
  ENDIF()
ENDIF()

ADD_LIBRARY( ${SPLINES} STATIC ${SOURCES} ${HEADERS} )
SET( SPLINES_LIB ${SPLINES} )

INSTALL( TARGETS ${TARGET}
         RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
         LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
         ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX} )

INSTALL( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include )

SET( BUILD_EXECUTABLE false )
IF( BUILD_EXECUTABLE )

  SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  FILE( MAKE_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin )

  MESSAGE( STATUS )
  MESSAGE( STATUS "GC_LIB = ${GC_LIB}" )
  MESSAGE( STATUS )

  FILE( GLOB S ./tests/*.cc )
  FOREACH (F ${S})
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    get_filename_component( TEXE ${RF} NAME_WE ) 
    ADD_EXECUTABLE( ${TEXE} ${RF} )
    TARGET_LINK_LIBRARIES( ${TEXE} ${SPLINES_LIB} ${GC_LIB} )
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
MESSAGE( STATUS "Sources = ${SOURCES}" )
MESSAGE( STATUS "Headers = ${HEADERS}" )
MESSAGE( STATUS "SPLINES = ${SPLINES}" )
MESSAGE( STATUS "GC_LIB  = ${GC_LIB}" )
