############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Università degli Studi di Trento                                    #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

project( Splines )

if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/cmake_utils")
  set( HOME "${PROJECT_SOURCE_DIR}/cmake_utils/" )
else()
  if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/../cmake_utils")
    set( HOME "${PROJECT_SOURCE_DIR}/../cmake_utils/" )
  else()
    set( HOME "${PROJECT_SOURCE_DIR}/../../cmake_utils/" )
  endif()
endif()


include( "${HOME}/CMakeLists-common.txt" )

project(
  ${PROJECT_NAME}
  VERSION ${UTILS_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
)

if (APPLE)
  #set( CMAKE_OSX_DEPLOYMENT_TARGET "13.3" )
  execute_process( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
  set( CMAKE_OSX_ARCHITECTURES ${ARCHITECTURE} )
endif()

include( "${HOME}/CMakeLists-cflags.txt" )
include( "${HOME}/CMakeLists-utilities.txt" )

# Evaluating the suffix to append to target built by this CMake file
utils_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  utils_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

if( EXISTS "${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt" )
  include(${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#
set( SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )

set( SOURCES )
file( GLOB S ${SOURCE_DIR}/*.cc ${SOURCE_DIR}/**/*.cc )
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})
list(FILTER SOURCES EXCLUDE REGEX ".*hsl_fake.cc$")

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB S
  ${SOURCE_DIR}/*.h*
  ${SOURCE_DIR}/*/*.h*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

set( PRIVATE_HEADERS )
# no private headers

list(FILTER PUBLIC_HEADERS EXCLUDE REGEX ".*.tmpl$")
list(FILTER PRIVATE_HEADERS EXCLUDE REGEX ".*.tmpl$")

# extra include directory for compilation
include_directories( ./lib3rd/include )

#       _                           _                 _
#    __| | ___ _ __   ___ _ __   __| | ___ _ __   ___(_) ___  ___
#   / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| |/ _ \/ __|
#  | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |  __/\__ \
#   \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|_|\___||___/
#             |_|
#
set(DEPEND_TARGETS)

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
utils_add_library_or_submodule( "quarticRootsFlocke" ROOTS )
utils_add_library_or_submodule( "UtilsLite" UTILS )
utils_add_library_or_submodule( "GenericContainer" GC "GC" )

message( STATUS "----------------------------------------------------------------" )
message( STATUS "UTILS = ${UTILS}" )
message( STATUS "ROOTS = ${ROOTS}" )
message( STATUS "GC    = ${GC}" )
message( STATUS "----------------------------------------------------------------" )

add_definitions( -DSPLINES_NO_COMPATIBILITY )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
utils_setup_target(
  ${PROJECT_NAME}
  TRUE # build static
  "${SOURCE_DIR}"
  "${SOURCES}"
  "${SOURCE_DIR}"
  "${PUBLIC_HEADERS}"
)
if ( UTILS_BUILD_SHARED )
  utils_setup_target(
    ${PROJECT_NAME}
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "${SOURCES}"
    "${SOURCE_DIR}"
    "${PUBLIC_HEADERS}"
  )
  target_link_libraries( ${UTILS_NAMESPACE}_${PROJECT_NAME} ${ROOTS} ${GC} ${UTILS} )
endif()

if ( DEPEND_TARGETS )
  add_dependencies( ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${DEPEND_TARGETS} )
endif()

#   _____       _
#  |_   _|__ __| |_ ___
#    | |/ -_|_-<  _(_-<
#    |_|\___/__/\__/__/
#
if ( UTILS_ENABLE_TESTS )

  enable_testing()

  set( EXELISTBIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin" )
  set( EXELISTCPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src_tests" )

  set(
    EXELISTCPP
    test1 test2 test3 test4 test5 test6 test8 test9 test10
  )

  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  set( LIBS ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${UTILS} ${ROOTS} ${GC} ${CMAKE_DL_LIBS} )
  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc )
    target_link_libraries( ${S} ${LIBS} )
    set_target_properties( ${S} PROPERTIES SUFFIX ".exe" )
    add_test( NAME "${S}" COMMAND ./bin/${S}.exe WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
  endforeach()

endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  # EXPORT (da definire)
  LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  PUBLIC_HEADER DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
)

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

# do not work
#include( ./cmake_utils/CMakeLists-cpack.txt )

utils_final_messages()
