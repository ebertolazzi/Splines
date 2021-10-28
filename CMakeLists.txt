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

cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_VERBOSE_MAKEFILE OFF)
# set(CMAKE_INSTALL_MESSAGE NEVER)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

if(ENABLE_CTEST)
  enable_testing()
endif()

# get version from the git tag
execute_process(COMMAND git describe --tags --abbrev=0 OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "\n$" "" VERSION "${VERSION}")


project(
  Splines
  VERSION ${VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
)

# NOTE: Predictable target name
set(STATIC_PROJECT "${PROJECT_NAME}Static")
set(DYNAMIC_PROJECT "${PROJECT_NAME}")

include(./CPackConfig.cmake)
include(./CMakeLists-cflags.txt)
include(./CMakeUtilities.txt)

# Evaluating the suffix to append to target built by this CMake file
ebertolazzi_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if(BUILD_SHARED)
  ebertolazzi_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
set(UPDIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(DIRGC    "${CMAKE_CURRENT_SOURCE_DIR}/extern/GC")
set(DIRUTILS "${CMAKE_CURRENT_SOURCE_DIR}/extern/Utils")
set(DIRROOTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/Flocke")
if(FALSE AND EXISTS "${UPDIR}/GC")
  add_subdirectory("${UPDIR}/GC" "${DIRGC}")
else()
  #add_subdirectory(./submodules/GenericContainer "${DIRGC}") # NON FUNZIONA
  add_subdirectory(./submodules/GenericContainer)
endif()
if(FALSE AND EXISTS "${UPDIR}/Utils")
  add_subdirectory("${UPDIR}/Utils" "${DIRUTILS}")
else()
  #add_subdirectory(./submodules/Utils "${DIRUTILS}") # NON FUNZIONA
  add_subdirectory(./submodules/Utils)
endif()
if(FALSE AND EXISTS "${UPDIR}/quarticRootsFlocke")
  add_subdirectory("${UPDIR}/quarticRootsFlocke" "${DIRROOTS}")
else()
  #add_subdirectory(./submodules/quarticRootsFlocke "${DIRROOTS}") # NON FUNZIONA
  add_subdirectory(./submodules/quarticRootsFlocke)
endif()

include_directories(lib3rd/include)

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
set(CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG})

# NOTE: Glob shall be avoided. Headers and sources shall be declared 
# manually in order to rebuild cmake-cache automatically.
# Ref: https://cmake.org/cmake/help/v3.15/command/file.html?highlight=glob#filesystem
set(Splines_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

set( SOURCES )
file( GLOB S ./src/*.cc )
foreach (F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach (F ${S})

set( HEADERS )
FILE( GLOB S ./src/*.hh )
foreach (F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND HEADERS ${RF} )
endforeach (F ${S})

FILE( GLOB S ./src/*.h* )
foreach (F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND HEADERS ${RF} )
endforeach (F ${S})

set( HEADERS2 )
FILE( GLOB S ./src/Splines/*.h* )
foreach (F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND HEADERS2 ${RF} )
endforeach (F ${S})


add_library(${STATIC_PROJECT} STATIC ${SOURCES})
target_include_directories(
  ${STATIC_PROJECT} PUBLIC 
  "$<BUILD_INTERFACE:${Splines_DIR}>"
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>"
)
set_target_properties(
  ${STATIC_PROJECT} PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}${ARTIFACTS_STATIC_SUFFIX}"
  # FIXME: I'm not sure if this is required for OSX Framework
  #PUBLIC_HEADER "${HEADERS}"
  POSITION_INDEPENDENT_CODE ON
)

if(BUILD_SHARED)
  add_library(${DYNAMIC_PROJECT} SHARED ${SOURCES})
  set_target_properties(
    ${DYNAMIC_PROJECT} PROPERTIES
    OUTPUT_NAME "${PROJECT_NAME}${ARTIFACTS_DYNAMIC_SUFFIX}"
    # FIXME: I'm not sure if this is required for OSX Framework
    # PUBLIC_HEADER "${HEADERS}"
    POSITION_INDEPENDENT_CODE ON
  )
  target_include_directories(
    ${DYNAMIC_PROJECT} PUBLIC 
    "$<BUILD_INTERFACE:${Splines_DIR}>"
    "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>"
  )
  target_link_libraries( ${DYNAMIC_PROJECT} PUBLIC GenericContainerStatic UtilsStatic QuarticStatic )
  if(WIN32)
    set_target_properties(
      ${DYNAMIC_PROJECT} PROPERTIES 
      WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
  endif()
endif()

if(BUILD_EXECUTABLE OR ENABLE_CTEST)
  enable_testing()
  set(TEST_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")
  set(
    EXECUTABLES
    test1 test2 test3 test4 test5
    test6 test8 test9 test10
  )
  foreach(EXE ${EXECUTABLES})
    add_executable(${EXE} tests-cpp/${EXE}.cc)
    target_link_libraries(${EXE} PRIVATE ${STATIC_PROJECT} QuarticStatic UtilsStatic GenericContainerStatic ${CMAKE_DL_LIBS})
    set_target_properties(
      ${EXE} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${TEST_OUTPUT_PATH}"
    )
    add_test(${EXE} "${CMAKE_CURRENT_SOURCE_DIR}/bin/${EXE}")
  endforeach()
endif()

# INSTALL TARGETS

set(INSTALL_TARGETS ${STATIC_PROJECT})
add_dependencies( ${STATIC_PROJECT} QuarticStatic UtilsStatic GenericContainerStatic )
if(BUILD_SHARED)
  add_dependencies( ${STATIC_PROJECT} QuarticStatic UtilsStatic GenericContainerStatic )
  set(INSTALL_TARGETS ${INSTALL_TARGETS} ${DYNAMIC_PROJECT})
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeInstall.txt")

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
