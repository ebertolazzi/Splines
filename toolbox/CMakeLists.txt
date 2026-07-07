#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  version: 1.0   date 9/3/2022                                         #
#                                                                       #
#  Copyright (C) 2022                                                   #
#                                                                       #
#      Enrico Bertolazzi, Francesco Biral and Paolo Bosetti             #
#      Dipartimento di Ingegneria Industriale                           #
#      Università degli Studi di Trento                                 #
#      Via Sommarive 9, I-38123, Trento, Italy                          #
#      email: enrico.bertolazzi@unitn.it                                #
#             francesco.biral@unitn.it                                  #
#             paolo.bosetti@unitn.it                                    #
#-----------------------------------------------------------------------#

cmake_minimum_required( VERSION 3.14 )

option( UTILS_GIT_VERSION   "Get version from GIT command"           ON )
option( UTILS_ENABLE_TESTS  "Enable compilation of testing targets"  ON )
option( UTILS_BUILD_SHARED  "Enable shared option compilation"       ON )
option( UTILS_INSTALL_LOCAL "Enable installation on local directory" ON )
option( UTILS_TWEAK_VERSION_IS_TIMESTAMP "Sets tweak version to timestamp" OFF )

set( TARGETS_TO_BE_INSTALLED )
set( CMAKE_INSTALL_MESSAGE         NEVER )
set( CMAKE_EXPORT_COMPILE_COMMANDS TRUE )
set( CMAKE_VERBOSE_MAKEFILE        FALSE )
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_INSTALL_RPATH $ORIGIN )
set( CMAKE_BUILD_WITH_INSTALL_RPATH ON )

if( NOT UNIX )
  set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
endif()

add_compile_definitions(_USE_MATH_DEFINES)

set(
  UTILS_CURRENT_VERSION "1.0"
  CACHE STRING "Default project version, if not taken from git"
)

if ( UTILS_GIT_VERSION AND NOT WIN32 )
  execute_process(
    COMMAND git describe --tags --abbrev=0
    OUTPUT_VARIABLE UTILS_PROJECT_VERSION
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

if ( NOT UTILS_PROJECT_VERSION )
  set( UTILS_BUILD_TIMESTAMP "0" )
  if( UTILS_TWEAK_VERSION_IS_TIMESTAMP )
    string(TIMESTAMP UTILS_BUILD_TIMESTAMP "%y%m%d%H")
  endif()
  set( UTILS_PROJECT_VERSION "${UTILS_CURRENT_VERSION}.${UTILS_BUILD_TIMESTAMP}" )
endif()

project( Splines VERSION ${UTILS_PROJECT_VERSION} )

if ( UTILS_INSTALL_LOCAL )
  set( CMAKE_INSTALL_PREFIX ${CMAKE_HOME_DIRECTORY} )
endif()

message("${PROJECT_NAME} current version = v${UTILS_PROJECT_VERSION}")

if ( NOT DEFINED CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
endif()

if ( NOT Matlab_ROOT_DIR AND DEFINED ENV{Matlab_ROOT_DIR} AND NOT "$ENV{Matlab_ROOT_DIR}" STREQUAL "" )
  set( Matlab_ROOT_DIR "$ENV{Matlab_ROOT_DIR}" )
endif()

if ( APPLE AND NOT Matlab_ROOT_DIR )
  file( GLOB MATLAB_APP_BUNDLES "/Applications/MATLAB_R*.app" )
  if ( MATLAB_APP_BUNDLES )
    list( SORT MATLAB_APP_BUNDLES ORDER DESCENDING )
    list( GET MATLAB_APP_BUNDLES 0 Matlab_ROOT_DIR )
  endif()
endif()

if ( Matlab_ROOT_DIR )
  set( Matlab_ROOT_DIR "${Matlab_ROOT_DIR}" CACHE PATH "MATLAB root directory" FORCE )
  message( STATUS "Using Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
else()
  message( STATUS "Matlab_ROOT_DIR not preset, relying on CMake autodetection" )
endif()

find_package( Matlab REQUIRED COMPONENTS MX_LIBRARY MEX_LIBRARY MAIN_PROGRAM )

# Aggiunge le directory di include di MATLAB
include_directories(${Matlab_INCLUDE_DIRS})

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )
message( STATUS "Matlab_MEX_LIBRARY = ${Matlab_MEX_LIBRARY}" )
message( STATUS "Matlab_MX_LIBRARY = ${Matlab_MX_LIBRARY}" )

if ( (UNIX OR LINUX) AND NOT APPLE )
  find_package( Threads REQUIRED )
  set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined" )
  set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined" )
endif()

set(COMMON_CXX_FLAGS
  "-pedantic"
  "-Wall"
  "-Wextra"
  "-Wshadow"
  "-Wnon-virtual-dtor"
  "-Wold-style-cast"
  "-Wcast-align"
  "-Wunused"
  "-Woverloaded-virtual"
  "-Wpedantic"
  "-Wconversion"
  "-Wsign-conversion"
  "-Wdouble-promotion"
  "-Wformat=2"
  "-Wimplicit-fallthrough"
  "-Wreturn-local-addr"
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND COMMON_CXX_FLAGS
    "-Wmisleading-indentation"
    "-Wduplicated-cond"
    "-Wduplicated-branches"
    "-Wlogical-op"
    "-Wnull-dereference"
    "-Wuseless-cast"
  )
endif()

string(JOIN " " COMMON_CXX_FLAGS_STR ${COMMON_CXX_FLAGS})

set(MSVC_COMPILE_FLAGS
  "/permissive"
  "/W4"
  "/w14242"
  "/w14254"
  "/w14263"
  "/w14265"
  "/w14287"
  "/we4289"
  "/w14296"
  "/w14311"
  "/w14545"
  "/w14546"
  "/w14547"
  "/w14549"
  "/w14555"
  "/w14619"
  "/w14640"
  "/w14826"
  "/w14905"
  "/w14906"
  "/w14928"
  "/external:W0"
  "/external:anglebrackets"
  "/external:templates-compat"
)

string(JOIN " " MSVC_COMPILE_FLAGS_STR ${MSVC_COMPILE_FLAGS})

set( VS_DEFINE  "/D_WINDOWS /D_SCL_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC" )
set( VS_FLAGS   "/nologo /std:c++17 /utf-8 /GS /Gy /GF /EHsc /bigobj" )

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
  set( CMAKE_CXX_FLAGS         "-std=c++17 ${COMMON_CXX_FLAGS_STR}" )
  set( CMAKE_CXX_FLAGS_RELEASE "-fPIC -O2 -g " )
  set( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -O0 -ggdb " )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  set( CMAKE_CXX_FLAGS         "-std=c++17 -stdlib=libc++ ${COMMON_CXX_FLAGS_STR}" )
  set( CMAKE_CXX_FLAGS_RELEASE "-fPIC -Os -g -fcolor-diagnostics " )
  set( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -O0 -gfull -fcolor-diagnostics " )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
  set( CMAKE_CXX_FLAGS         " ${MSVC_COMPILE_FLAGS_STR} ${VS_DEFINE} /sdl" )
  set( CMAKE_CXX_FLAGS_RELEASE " ${VS_FLAGS} /O2 /MD /Zo " )
  set( CMAKE_CXX_FLAGS_DEBUG   " ${VS_FLAGS} /Od /Ob0 /MDd /guard:cf /Zi /DEBUG:FULL " )
else()
  message( FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}" )
endif()

set( CMAKE_C_FLAGS         "" )
set( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
set( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

function(toolbox_artifacts_suffix VAR_NAME BUILD_STATIC)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(suffix "_debug")
  else()
    set(suffix "")
  endif()
  if ( ${BUILD_STATIC} )
    set(suffix "_static${suffix}")
  endif()
  if(APPLE)
    set( ${VAR_NAME} _osx${suffix} PARENT_SCOPE )
  elseif(MINGW)
    set( ${VAR_NAME} _mingw${suffix} PARENT_SCOPE )
  elseif(UNIX)
    set( ${VAR_NAME} _linux${suffix} PARENT_SCOPE )
  else()
    set( ${VAR_NAME} "_win_${BITS}${suffix}" PARENT_SCOPE )
  endif()
endfunction()

function(toolbox_final_messages)
  message(STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}")
  message(STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}")
  message(STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}")
  message(STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}")
  message(STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}")
  message(STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}")
  message(STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}")
  message(STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}")
  if(APPLE)
    message( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
    message( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
    message( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
  endif()
  message(STATUS "SOURCES = ${SOURCES}")
  message(STATUS "HEADERS = ${HEADERS}")
  message(STATUS "TARGETS_TO_BE_INSTALLED  = ${TARGETS_TO_BE_INSTALLED}")
  message(STATUS "CMAKE_HOME_DIRECTORY     = ${CMAKE_HOME_DIRECTORY}")
  message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
  message(STATUS "CMAKE_INSTALL_PREFIX     = ${CMAKE_INSTALL_PREFIX}")
  message(STATUS "UTILS_GIT_VERSION   = ${UTILS_GIT_VERSION}")
  message(STATUS "UTILS_ENABLE_TESTS  = ${UTILS_ENABLE_TESTS}")
  message(STATUS "UTILS_BUILD_SHARED  = ${UTILS_BUILD_SHARED}")
  message(STATUS "UTILS_INSTALL_LOCAL = ${UTILS_INSTALL_LOCAL}")
endfunction()

# Evaluating the suffix to append to target built by this CMake file
toolbox_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  toolbox_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ____                                _     _ _
#  |  _ \ _ __ ___ _ __ ___  __ _ _   _(_)___(_) |_ ___
#  | |_) | '__/ _ \ '__/ _ \/ _` | | | | / __| | __/ _ \
#  |  __/| | |  __/ | |  __/ (_| | |_| | \__ \ | ||  __/
#  |_|   |_|  \___|_|  \___|\__, |\__,_|_|___/_|\__\___|
#                              |_|
# MUST be executed before the build of the header list

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#

set( SOURCES )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src_lib/*.cc
)
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})

list(FILTER SOURCES EXCLUDE REGEX "Utils_nonlinear_system_tests\\.cc$")

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h*
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PUBLIC_HEADERS EXCLUDE REGEX "\.DS")

set( PRECOMPILED_HEADERS )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hh
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}" )
  list( APPEND PRECOMPILED_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PRECOMPILED_HEADERS EXCLUDE REGEX "\.DS")

# extra include directory for compilation
include_directories( src src/Splines src_lib src/Utils src/Utils/3rd )

message( STATUS "SOURCES = ${SOURCES}" )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#

toolbox_final_messages()

add_library( Splines STATIC ${SOURCES} )

# GenericContainer now requires C++20.
set_target_properties(Splines PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

set(
  MEX_LIST
  BaseHermiteMexWrapper
  Spline1DMexWrapper
  Spline1DblendMexWrapper
  Spline2DMexWrapper
  Spline2DblendMexWrapper
  SplineSetMexWrapper
  SplineVecMexWrapper
)

message( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if ( CMAKE_BUILD_TYPE MATCHES Debug )
  message( STATUS "DEBUG MODE\n\n")
  add_compile_definitions(DEBUG)
endif()

if ( (UNIX OR LINUX) AND NOT APPLE )
  set( EXTRA_LIBS Threads::Threads ${CMAKE_DL_LIBS} )
else()
  set( EXTRA_LIBS ${CMAKE_DL_LIBS} )
endif()

set( EXTRA_LIBS ${EXTRA_LIBS} Splines )

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Set C++ standard to 20 for C++ files only
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# ----------------------------------------------------------------------
# Calcola automaticamente il suffisso MEX per la piattaforma corrente
# ----------------------------------------------------------------------
if(APPLE)
    # macOS: determina se è Apple Silicon (ARM64) o Intel
    execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE ARCHITECTURE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(ARCHITECTURE STREQUAL "arm64" OR ARCHITECTURE STREQUAL "aarch64")
        set(MEX_SUFFIX ".mexmaca64")
        message(STATUS "Detected Apple Silicon (ARM64), using MEX suffix: ${MEX_SUFFIX}")
    else()
        set(MEX_SUFFIX ".mexmaci64")
        message(STATUS "Detected macOS Intel, using MEX suffix: ${MEX_SUFFIX}")
    endif()
elseif(WIN32)
    set(MEX_SUFFIX ".mexw64")
    message(STATUS "Detected Windows, using MEX suffix: ${MEX_SUFFIX}")
elseif(UNIX AND NOT APPLE)
    set(MEX_SUFFIX ".mexa64")
    message(STATUS "Detected Linux, using MEX suffix: ${MEX_SUFFIX}")
else()
    set(MEX_SUFFIX ".mex")
    message(STATUS "Unknown platform, using default MEX suffix: ${MEX_SUFFIX}")
endif()

# ----------------------------------------------------------------------
# Funzione per creare target MEX con nome di output corretto
# ----------------------------------------------------------------------
function(add_mex_target MEX_NAME)
    message(STATUS "Building MEX: ${MEX_NAME}")
    
    # Crea il target MEX
    add_library(mex_${MEX_NAME} MODULE 
        src_mex/mex_${MEX_NAME}.cc
        ${Matlab_ROOT_DIR}/extern/version/c_mexapi_version.c
    )
    
    # Imposta le proprietà con OUTPUT_NAME che include il suffisso
    set_target_properties(mex_${MEX_NAME} PROPERTIES
        PREFIX ""
        OUTPUT_NAME "${MEX_NAME}${MEX_SUFFIX}"  # Nome completo con suffisso
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        BUNDLE TRUE
    )
    
    # Aggiungi le directory di include
    target_include_directories(mex_${MEX_NAME} 
        PRIVATE 
        ${Matlab_INCLUDE_DIRS}
        src 
        src/Splines
        src_lib 
        src/Utils
        src/Utils/3rd
    )
    
    # Aggiungi le definizioni di compilazione per MEX
    target_compile_definitions(mex_${MEX_NAME} PRIVATE
        MATLAB_MEX_FILE
        MATLAB_DEFAULT_RELEASE=R2018a
        _USE_MATH_DEFINES
    )
    
    # Link delle librerie
    target_link_libraries(mex_${MEX_NAME}
        ${EXTRA_LIBS}
        ${Matlab_MEX_LIBRARY}
        ${Matlab_MX_LIBRARY}
    )
    
    # Opzioni di linker specifiche per macOS
    if(APPLE)
        target_link_options(mex_${MEX_NAME} PRIVATE
            "-Wl,-undefined,dynamic_lookup"
            "-L${Matlab_ROOT_DIR}/bin/maca64"
            "-Wl,-rpath,${Matlab_ROOT_DIR}/bin/maca64"
        )
    endif()
    
    # Copia il file MEX nella directory bin
    add_custom_command(
        TARGET mex_${MEX_NAME}
        POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:mex_${MEX_NAME}>"
        "${CMAKE_CURRENT_SOURCE_DIR}/bin"
        COMMENT "Copying ${MEX_NAME}${MEX_SUFFIX} to bin directory"
    )
endfunction()

# ----------------------------------------------------------------------
# Crea tutti i target MEX
# ----------------------------------------------------------------------
foreach(MEX ${MEX_LIST})
    add_mex_target(${MEX})
endforeach()

# Aggiungere un target per la pulizia
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMENT "Cleaning all build files"
)

# Aggiungere informazioni di debug
message(STATUS "Matlab libraries found:")
message(STATUS "  MEX Library: ${Matlab_MEX_LIBRARY}")
message(STATUS "  MX Library: ${Matlab_MX_LIBRARY}")
message(STATUS "  Include dirs: ${Matlab_INCLUDE_DIRS}")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "  MEX suffix: ${MEX_SUFFIX}")
