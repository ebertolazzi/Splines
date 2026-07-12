############################################################################
#                                                                          #
#  Splines -- modern, idiomatic CMake build.                              #
#                                                                          #
#  Self-contained: no dependency on the external `cmake_utils` grab-bag   #
#  of .cmake fragments (it isn't even a proper CMake subproject -- it has #
#  no CMakeLists.txt of its own -- and its Clang flags hardcode           #
#  `-stdlib=libc++`, which breaks on Linux unless libc++-dev happens to   #
#  be installed). Structure mirrors GenericContainer's CMakeLists.txt,    #
#  which this project itself depends on.                                 #
#                                                                          #
############################################################################

cmake_minimum_required( VERSION 3.25 )

# --------------------------------------------------------------------------
# Apple FAT/universal build setup
# --------------------------------------------------------------------------
# On Apple platforms a single library can contain slices for several CPU
# architectures, for example arm64 and x86_64.  This is controlled by
# CMAKE_OSX_ARCHITECTURES.
#
# The variable must be set before project(), otherwise CMake may initialize
# compilers and targets before the requested architectures are known.
#
# Default behavior:
#   - ON  on Apple hosts or Apple cross-compilation targets.
#   - OFF on Linux/Windows and other platforms, where there is no equivalent
#     single-file universal library format.
#
# User overrides:
#   -DUTILS_BUILD_FAT=OFF
#   -DUTILS_FAT_ARCHITECTURES=arm64
#   -DCMAKE_OSX_ARCHITECTURES=arm64\;x86_64
set(UTILS_BUILD_FAT_DEFAULT OFF)
if(CMAKE_HOST_APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin|iOS|tvOS|watchOS|visionOS")
  set(UTILS_BUILD_FAT_DEFAULT ON)
endif()
option(UTILS_BUILD_FAT "Build Apple universal/fat binaries when supported" ${UTILS_BUILD_FAT_DEFAULT})
set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")

if(UTILS_BUILD_FAT AND NOT CMAKE_OSX_ARCHITECTURES)
  set(
    CMAKE_OSX_ARCHITECTURES
    "${UTILS_FAT_ARCHITECTURES}"
    CACHE STRING "Build architectures for Apple platforms"
    FORCE
  )
endif()


# --------------------------------------------------------------------------
# Project version discovery
# --------------------------------------------------------------------------
# Prefer the most recent Git tag as the project version when the source tree is
# a Git checkout.  This keeps the CMake package version aligned with release
# tags without requiring a manually edited version string in this file.
#
# Tags that start with a leading "v" are normalized by removing that prefix,
# so both "v1.2.3" and "1.2.3" become a valid CMake project version.  When Git
# is not installed, or the directory is not a tagged Git checkout, the build
# falls back to "0.0.0" so configuration remains deterministic.
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(GIT_VERSION MATCHES "^v")
    string(SUBSTRING ${GIT_VERSION} 1 -1 GIT_VERSION)
  endif()
endif()
if(NOT GIT_VERSION)
  set(GIT_VERSION "0.0.0")
endif()

project(
  Splines
  VERSION ${GIT_VERSION}
  DESCRIPTION "A C++ library for univariate and bivariate spline interpolation"
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
  LANGUAGES CXX C
)

# ============================================================================
# Platform suffix used in generated library filenames
# ============================================================================
# The real library artifacts receive an explicit platform suffix, e.g.
#   libSplines_osx_static.a
#   libSplines_linux_static.a
#   Splines_win64.dll / import library variants on Windows
#
# A separate compatibility alias with the generic name is created later.  The
# suffix makes artifacts from different platforms distinguishable when copied
# into the same local lib directory or distributed as prebuilt binaries.
if ( APPLE )
  set( SPLINES_PLATFORM_SUFFIX osx )
elseif ( WIN32 )
  if ( MINGW )
    if ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
      set( SPLINES_PLATFORM_SUFFIX mingw64 )
    else()
      set( SPLINES_PLATFORM_SUFFIX mingw32 )
    endif()
  elseif ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( SPLINES_PLATFORM_SUFFIX win64 )
  else()
    set( SPLINES_PLATFORM_SUFFIX win32 )
  endif()
elseif ( UNIX )
  set( SPLINES_PLATFORM_SUFFIX linux )
else()
  string( TOLOWER "${CMAKE_SYSTEM_NAME}" SPLINES_PLATFORM_SUFFIX )
endif()

# ============================================================================
# User-facing options and basic build defaults
# ============================================================================
# This section defines the configuration switches intended for users, package
# scripts, CI jobs, and parent projects using FetchContent.  Defaults are chosen
# so the project behaves as a complete standalone build at the top level, while
# avoiding tests/install side effects when consumed as a dependency.

# Default to a Release build when the generator uses a single configuration
# and the user hasn't picked one (multi-config generators like Ninja
# Multi-Config or Visual Studio ignore CMAKE_BUILD_TYPE entirely).
if ( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES )
  set( CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE )
  set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel" )
endif()

# Standalone vs. consumed-as-a-dependency: tests only make sense standalone,
# and installing shouldn't be forced on a project that FetchContent's us.
option( BUILD_SHARED_LIBS "Build Splines as a shared library" OFF )
option( BUILD_TESTING "Build the Splines test executables" ${PROJECT_IS_TOP_LEVEL} )
option( SPLINES_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL} )
option( SPLINES_BUILD_BENCHMARKS "Build the micro-benchmarks" ${PROJECT_IS_TOP_LEVEL} )
option( SPLINES_STRICT_WARNINGS "Compile the Splines target with narrowing/conversion warnings" OFF )
option( SPLINES_UPDATE_3RDPARTY "Refresh lib3rd/include from resolved dependencies at configure time" ${PROJECT_IS_TOP_LEVEL} )
option( SPLINES_COLLECT_DEPENDENCIES "Deprecated alias for SPLINES_UPDATE_3RDPARTY" ${SPLINES_UPDATE_3RDPARTY} )
if ( SPLINES_COLLECT_DEPENDENCIES )
  set(
    SPLINES_UPDATE_3RDPARTY ON
    CACHE BOOL "Refresh lib3rd/include from resolved dependencies at configure time"
    FORCE
  )
endif()

# Local sibling checkout paths.  When these directories exist, they are used
# directly; otherwise the dependency is downloaded with FetchContent.
#
# This layout supports the usual development workspace where Splines,
# GenericContainer, UtilsLite, quarticRootsFlocke, Eigen, and nlohmann/json live
# next to each other.  Keeping these as CACHE paths lets users override them
# from the command line without editing this file.
set( SPLINES_GC_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../GC" CACHE PATH "Local GenericContainer checkout" )
set( SPLINES_UTILSLITE_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite" CACHE PATH "Local UtilsLite checkout" )
set( SPLINES_QUARTIC_ROOTS_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../quarticRootsFlocke" CACHE PATH "Local quarticRootsFlocke checkout" )
set( SPLINES_EIGEN_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Eigen" CACHE PATH "Local Eigen checkout" )
set( SPLINES_JSON_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../json" CACHE PATH "Local nlohmann_json checkout" )
set( SPLINES_EIGEN_VERSION "5.0.0" CACHE STRING "Eigen version to fetch" )

# Default install layout requested for this project:
#   ${CMAKE_CURRENT_SOURCE_DIR}/lib/include
#   ${CMAKE_CURRENT_SOURCE_DIR}/lib/lib
#
# The historical build layout writes headers and libraries under the source
# tree's lib/ directory.  Keeping that behavior here avoids breaking existing
# Rake/build scripts and downstream projects that expect this local layout.
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE )
endif()
set( CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE )
set( CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE )

# Do not set global CMAKE_*_OUTPUT_DIRECTORY variables here.  Output locations
# are assigned per target below, which avoids unintentionally redirecting
# artifacts produced by FetchContent dependencies.

include( GNUInstallDirs )
include( FetchContent )
include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Update3rdParties.cmake" )

# ============================================================================
# Dependency helper: prefer a local checkout, otherwise use FetchContent
# ============================================================================
# splines_use_local_or_fetch(NAME LOCAL_DIR) centralizes the policy used for
# third-party dependencies:
#   1. If LOCAL_DIR contains a CMake project, force FetchContent to use it.
#   2. Otherwise, let FetchContent populate the dependency from the declaration.
#   3. Record *_SOURCE and *_SOURCE_DIR variables for the final summary and for
#      optional third-party header collection.
#
# The helper temporarily disables install rules while adding dependencies.  This
# prevents dependency install() commands from being injected into Splines'
# install target and keeps this project responsible for its own install layout.
function( splines_use_local_or_fetch NAME LOCAL_DIR )
  # Normalize the user-provided local dependency path relative to this project.
  # This allows cache values such as ../UtilsLite to work independently of the
  # directory from which CMake was invoked.
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  string( TOUPPER "${NAME}" _name_upper )
  string( TOLOWER "${NAME}" _name_lower )
  # Reuse FetchContent's base directory when it is already configured; otherwise
  # follow CMake's conventional _deps location under the build tree.
  if ( FETCHCONTENT_BASE_DIR )
    set( _fetchcontent_base_dir "${FETCHCONTENT_BASE_DIR}" )
  else()
    set( _fetchcontent_base_dir "${CMAKE_BINARY_DIR}/_deps" )
  endif()
  set( _dependency_install_stub "${_fetchcontent_base_dir}/${_name_lower}-build/cmake_install.cmake" )
  set( _saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}" )
  set( CMAKE_SKIP_INSTALL_RULES ON )

  # A local checkout wins over downloading.  Setting FETCHCONTENT_SOURCE_DIR_* in
  # the cache is the documented FetchContent override mechanism.
  if ( EXISTS "${_local_abs}/CMakeLists.txt" )
    message( STATUS "Using local ${NAME}: ${_local_abs}" )
    set( "FETCHCONTENT_SOURCE_DIR_${_name_upper}" "${_local_abs}" CACHE PATH "Local ${NAME} checkout" FORCE )
    FetchContent_MakeAvailable( ${NAME} )
    set( "${NAME}_SOURCE" "local: ${_local_abs}" PARENT_SCOPE )
    set( "${NAME}_SOURCE_DIR" "${_local_abs}" PARENT_SCOPE )
  else()
    message( STATUS "Fetching ${NAME}" )
    FetchContent_MakeAvailable( ${NAME} )
    set( "${NAME}_SOURCE" "FetchContent" PARENT_SCOPE )
    string( TOLOWER "${NAME}" _name_lower )
    if ( DEFINED "${_name_lower}_SOURCE_DIR" )
      set( "${NAME}_SOURCE_DIR" "${${_name_lower}_SOURCE_DIR}" PARENT_SCOPE )
    endif()
  endif()

  # Some dependency projects expect a generated cmake_install.cmake in their
  # build directory.  Create a harmless stub because their install rules are
  # intentionally disabled while they are embedded in Splines.
  file( MAKE_DIRECTORY "${_fetchcontent_base_dir}/${_name_lower}-build" )
  if ( NOT EXISTS "${_dependency_install_stub}" )
    file(
      WRITE "${_dependency_install_stub}"
      "# Install disabled for dependency ${NAME} inside Splines.\n"
    )
  endif()

  set( CMAKE_SKIP_INSTALL_RULES "${_saved_skip_install_rules}" )
endfunction()

#   ___                          _                 _
#  |   \ ___ _ __  ___ _ _  __| |___ _ _  __ ___(_)___ ___
#  | |) / -_) '_ \/ -_) ' \/ _` / -_) ' \/ _/ -_)_(_-</ -_)
#  |___/\___| .__/\___|_||_\__,_\___|_||_\__\___(_)__/\___|
#            |_|
#
# nlohmann/json dependency.
#
# Splines links against the canonical nlohmann_json::nlohmann_json target.
# The local-or-FetchContent helper lets a developer use a sibling checkout
# while CI or clean builds can download the pinned release automatically.
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( nlohmann_json "${SPLINES_JSON_LOCAL_DIR}" )
set( SPLINES_NLOHMANN_JSON_SOURCE "${nlohmann_json_SOURCE}" )



# ============================================================================
# Eigen dependency
# ============================================================================
# Eigen is header-only for this project.  Resolution order is:
#   1. Reuse an existing Eigen3::Eigen target supplied by a parent project.
#   2. Use a local checkout from SPLINES_EIGEN_LOCAL_DIR.
#   3. Download the requested Eigen release with FetchContent.
#
# The resolved include root is stored in SPLINES_EIGEN_SOURCE_DIR because the
# optional lib3rd/include refresh step needs to copy Eigen headers.
set( SPLINES_EIGEN_SOURCE "FetchContent" )

if ( NOT TARGET Eigen3::Eigen )
  get_filename_component(
    SPLINES_EIGEN_LOCAL_DIR_ABS
    "${SPLINES_EIGEN_LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  if ( EXISTS "${SPLINES_EIGEN_LOCAL_DIR_ABS}/Eigen/Core" )
    message( STATUS "Using local Eigen checkout: ${SPLINES_EIGEN_LOCAL_DIR_ABS}" )
    set( SPLINES_EIGEN_SOURCE "local: ${SPLINES_EIGEN_LOCAL_DIR_ABS}" )
    set( SPLINES_EIGEN_SOURCE_DIR "${SPLINES_EIGEN_LOCAL_DIR_ABS}" )
    add_library( Eigen3::Eigen INTERFACE IMPORTED GLOBAL )
    target_include_directories( Eigen3::Eigen SYSTEM INTERFACE "${SPLINES_EIGEN_LOCAL_DIR_ABS}" )
  else()
    message( STATUS "Fetching Eigen ${SPLINES_EIGEN_VERSION}" )
    set( SPLINES_EIGEN_SOURCE "FetchContent: ${SPLINES_EIGEN_VERSION}" )
    FetchContent_Declare(
      eigen3
      URL "https://gitlab.com/libeigen/eigen/-/archive/${SPLINES_EIGEN_VERSION}/eigen-${SPLINES_EIGEN_VERSION}.tar.gz"
      SOURCE_SUBDIR cmake
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable( eigen3 )
    set( SPLINES_EIGEN_SOURCE_DIR "${eigen3_SOURCE_DIR}" )

    if ( NOT TARGET Eigen3::Eigen )
      add_library( Eigen3::Eigen INTERFACE IMPORTED GLOBAL )
      target_include_directories( Eigen3::Eigen SYSTEM INTERFACE "${eigen3_SOURCE_DIR}" )
    endif()
  endif()
else()
  message( STATUS "Using existing Eigen3::Eigen target" )
  set( SPLINES_EIGEN_SOURCE "existing target: Eigen3::Eigen" )
  get_target_property( _eigen_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES )
  if ( _eigen_inc )
    list( GET _eigen_inc 0 SPLINES_EIGEN_SOURCE_DIR )
  endif()
endif()

# ============================================================================
# GenericContainer dependency
# ============================================================================
# Splines uses the YAML-enabled GenericContainer target.  Install rules are
# disabled before bringing the dependency into the build so that a standalone
# Splines installation does not accidentally install GenericContainer as a
# separate package component.
set( GENERIC_CONTAINER_INSTALL OFF CACHE BOOL "Disable GenericContainer install rules inside Splines" FORCE )
FetchContent_Declare(
  GenericContainer
  GIT_REPOSITORY https://github.com/ebertolazzi/GenericContainer.git
  GIT_TAG        GC2.0
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( GenericContainer "${SPLINES_GC_LOCAL_DIR}" )
set( GC GenericContainer::Yaml )

# Clang standard-library propagation for GenericContainer.
#
# Some dependency targets may be INTERFACE libraries while others are concrete
# build targets.  Apply -stdlib=libc++ with the correct visibility in both
# cases so downstream compilation and linking stay consistent under Clang.
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  foreach( _gc_target GenericContainer GenericContainerYaml )
    if ( TARGET ${_gc_target} )
      get_target_property( _gc_type ${_gc_target} TYPE )
      if ( _gc_type STREQUAL "INTERFACE_LIBRARY" )
        target_compile_options( ${_gc_target} INTERFACE -stdlib=libc++ )
        target_link_options( ${_gc_target} INTERFACE -stdlib=libc++ )
      else()
        target_compile_options( ${_gc_target} PUBLIC -stdlib=libc++ )
        target_link_options( ${_gc_target} PUBLIC -stdlib=libc++ )
      endif()
    endif()
  endforeach()
endif()

# UtilsLite and quarticRootsFlocke share the same legacy option names for
# install rules and internal test suites.  Disable those features before adding
# the dependencies so this project remains the only owner of the install target
# and ctest only runs Splines' own tests.
set( UTILS_INSTALL_LOCAL OFF CACHE BOOL "Disable dependency-local install rules inside Splines" FORCE )
set( UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable dependency test suites" FORCE )

# ============================================================================
# UtilsLite dependency
# ============================================================================
# UtilsLite provides utility headers and library code used by Splines and its
# transitive dependencies.  It is added via the same local-or-FetchContent
# helper and then exposed through the canonical utils::UtilsLite_Static target.
FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        Utils2.0
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( UtilsLite "${SPLINES_UTILSLITE_LOCAL_DIR}" )
set( UTILS utils::UtilsLite_Static )

# Treat UtilsLite's headers as SYSTEM includes for anything that links it.
#
# This suppresses warnings from vendored or external dependency headers while
# preserving warnings for Splines' own sources.
foreach( _utils_tgt utils_UtilsLite_Static utils_UtilsLite )
  if ( TARGET ${_utils_tgt} )
    get_target_property( _utils_inc ${_utils_tgt} INTERFACE_INCLUDE_DIRECTORIES )
    if ( _utils_inc )
      set_target_properties(
        ${_utils_tgt} PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_utils_inc}"
      )
    endif()
  endif()
endforeach()

# ============================================================================
# quarticRootsFlocke dependency
# ============================================================================
# The spline algorithms depend on quartic root solving routines.  This target is
# pulled in as a normal build dependency and later moved to lib3rd/lib so the
# local artifact layout separates project outputs from dependency outputs.
FetchContent_Declare(
  quarticRootsFlocke
  GIT_REPOSITORY https://github.com/ebertolazzi/quarticRootsFlocke.git
  GIT_TAG        1.2.0
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( quarticRootsFlocke "${SPLINES_QUARTIC_ROOTS_LOCAL_DIR}" )
set( ROOTS utils_quarticRootsFlocke_Static )

set( SPLINES_DEPENDENCIES Eigen3::Eigen )

if ( SPLINES_UPDATE_3RDPARTY )
  splines_update_3rdparties(
    DESTINATION            "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include"
    EIGEN_ROOT             "${SPLINES_EIGEN_SOURCE_DIR}"
    JSON_ROOT              "${nlohmann_json_SOURCE_DIR}"
    GENERIC_CONTAINER_ROOT "${GenericContainer_SOURCE_DIR}"
    UTILSLITE_ROOT         "${UtilsLite_SOURCE_DIR}"
    QUARTIC_ROOTS_ROOT     "${quarticRootsFlocke_SOURCE_DIR}"
  )
endif()

message( STATUS "----------------------------------------------------------------" )
message( STATUS "UTILS = ${UTILS}" )
message( STATUS "ROOTS = ${ROOTS}" )
message( STATUS "GC    = ${GC} Eigen3::Eigen" )
message( STATUS "----------------------------------------------------------------" )

# ============================================================================
# Output directory helper for project targets
# ============================================================================
# Assigns archive, shared-library, and runtime output directories on a single
# target.  This is deliberately target-scoped: it avoids changing global CMake
# output variables and therefore avoids surprising dependencies added with
# FetchContent.
function( splines_set_target_output_dirs TARGET_NAME LIB_DIR )
  set_target_properties( ${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

# ============================================================================
# Output directory helper for dependency targets
# ============================================================================
# Moves concrete dependency artifacts under lib3rd/lib.  Alias targets are first
# resolved to their real target.  INTERFACE libraries are skipped because they do
# not produce binary artifacts.
function( splines_set_dependency_output_dirs LIB_DIR )
  # ARGN contains a variable-length list of dependency targets.  Each target may
  # be a real build target or an alias target, so resolve ALIASED_TARGET before
  # assigning output directories.
  foreach( tgt IN LISTS ARGN )
    get_target_property( real_tgt ${tgt} ALIASED_TARGET )
    if( NOT real_tgt )
      set( real_tgt ${tgt} )
    endif()
    if( TARGET ${real_tgt} )
      get_target_property( tgt_type ${real_tgt} TYPE )
      if( NOT tgt_type STREQUAL "INTERFACE_LIBRARY" )
        set_target_properties( ${real_tgt} PROPERTIES
          ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
          LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
        )
      endif()
    endif()
  endforeach()
endfunction()

#   _     _ _
#  | |   (_) |__  _ __ __ _ _ __ _   _
#  | |   | | '_ \| '__/ _` | '__| | | |
#  | |___| | |_) | | | (_| | |  | |_| |
#  |_____|_|_.__/|_|  \__,_|_|   \__, |
#                                |___/
#
# Main Splines library target.
#
# The source list is explicit rather than globbed so accidental files do not
# silently become part of the ABI.  The public alias Splines::Splines is the
# target that downstream CMake projects should link against.
add_library(
  Splines
  src/SplineAkima.cc
  src/SplineBiCubic.cc
  src/SplineBiCubicBase.cc
  src/SplineBilinear.cc
  src/SplineBiQuintic.cc
  src/SplineBiQuinticBase.cc
  src/SplineBuild.cc
  src/SplineConstant.cc
  src/SplineCubic.cc
  src/SplineCubicBase.cc
  src/SplineHermite.cc
  src/SplineLinear.cc
  src/SplineParametrization.cc
  src/SplinePchip.cc
  src/SplineQuintic.cc
  src/SplineQuinticBase.cc
  src/Splines.cc
  src/Splines1D.cc
  src/SplinesCinterface.cc
  src/SplineSet.cc
  src/SplineSurf.cc
  src/SplineUtils.cc
  src/SplineVanLeer.cc
  src/SplineVec.cc
)
add_library( Splines::Splines ALIAS Splines )

# Place the main project artifact in lib/lib, matching the repository's local
# build/install convention.
splines_set_target_output_dirs( Splines lib/lib )

target_compile_features( Splines PUBLIC cxx_std_20 )

# Public compile definitions become part of the usage requirements for clients
# linking Splines::Splines.  SPLINES_NO_COMPATIBILITY selects the modern API
# mode; _USE_MATH_DEFINES enables mathematical constants on MSVC-compatible
# standard library implementations.
target_compile_definitions( Splines PUBLIC SPLINES_NO_COMPATIBILITY )
target_compile_definitions( Splines PUBLIC _USE_MATH_DEFINES )
# The build interface points to the repository include directory, while the
# install interface is relative to the final install prefix.  This keeps the
# target relocatable for installed packages.
target_include_directories(
  Splines PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
# Link dependencies as PUBLIC because Splines headers expose types or compile
# definitions that require downstream consumers to see the same dependency
# usage requirements.
target_link_libraries( Splines PUBLIC ${ROOTS} ${GC} ${UTILS} nlohmann_json::nlohmann_json )

set_target_properties(
  Splines PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Splines
)

# Baseline warning policy for the project target.  Strict warning groups remain
# opt-in because they may be too noisy for some compiler/library combinations.
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
  target_compile_options( Splines PRIVATE -Wall -Wextra )
  if ( SPLINES_STRICT_WARNINGS )
    target_compile_options(
      Splines PRIVATE
      -Wconversion -Wsign-conversion -Wshadow -Wdouble-promotion
    )
  endif()
endif()

if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  target_compile_options( Splines PUBLIC -stdlib=libc++ )
  target_link_options( Splines PUBLIC -stdlib=libc++ )
endif()

# MSVC-specific options: /utf-8 avoids source encoding ambiguity, while /bigobj
# raises the object section limit for template-heavy translation units.
if ( MSVC )
  target_compile_options( Splines PUBLIC /utf-8 /bigobj )
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )

# ============================================================================
# Dependency artifact placement
# ============================================================================
# Dependencies that produce binary libraries are redirected to lib3rd/lib.  This
# keeps lib/lib focused on the Splines public artifact while still keeping all
# build products inside the source-tree local layout used by the project tools.
splines_set_dependency_output_dirs( lib3rd/lib
  GenericContainer
  Lua
  ${UTILS}
  ${ROOTS}
)

# ============================================================================
# Platform-suffixed artifact naming and compatibility aliases
# ============================================================================
# The real output is named with the platform suffix, for example:
#   libSplines_osx_static.a
#   libSplines_linux_static.a
#
# A generic compatibility name, such as libSplines.a, is then created as a
# symlink on Unix-like platforms or a copied file on Windows.  This preserves
# compatibility with older scripts that link by the generic name while still
# keeping the real file name platform-specific.
function( splines_add_platform_library_alias TARGET_NAME BASE_NAME )
  # Compute the real platform-specific filename and the generic compatibility
  # filename.  Shared-library builds need to account for runtime files and, on
  # Windows, import libraries; static builds only need the archive name.
  if ( BUILD_SHARED_LIBS )
    set( TARGET_OUTPUT_NAME "${BASE_NAME}_${SPLINES_PLATFORM_SUFFIX}" )

    if ( WIN32 )
      set( TARGET_REAL_RUNTIME "${CMAKE_SHARED_LIBRARY_PREFIX}${TARGET_OUTPUT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
      set( TARGET_LINK_RUNTIME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
      set( TARGET_REAL_IMPORT  "${CMAKE_IMPORT_LIBRARY_PREFIX}${TARGET_OUTPUT_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}" )
      set( TARGET_LINK_IMPORT  "${CMAKE_IMPORT_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_IMPORT_LIBRARY_SUFFIX}" )
    else()
      set( TARGET_REAL_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${TARGET_OUTPUT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
      set( TARGET_LINK_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
    endif()
  else()
    set( TARGET_OUTPUT_NAME "${BASE_NAME}_${SPLINES_PLATFORM_SUFFIX}_static" )

    if ( WIN32 )
      set( TARGET_REAL_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${TARGET_OUTPUT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
      set( TARGET_LINK_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    else()
      set( TARGET_REAL_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${TARGET_OUTPUT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
      set( TARGET_LINK_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    endif()
  endif()

  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    OUTPUT_NAME         "${TARGET_OUTPUT_NAME}"
    ARCHIVE_OUTPUT_NAME "${TARGET_OUTPUT_NAME}"
    LIBRARY_OUTPUT_NAME "${TARGET_OUTPUT_NAME}"
    RUNTIME_OUTPUT_NAME "${TARGET_OUTPUT_NAME}"
  )

  if ( WIN32 )
    if ( BUILD_SHARED_LIBS )
      add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:${TARGET_NAME}>"
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_RUNTIME}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_LINKER_FILE:${TARGET_NAME}>"
                "$<TARGET_LINKER_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_IMPORT}"
        VERBATIM
      )

      if ( SPLINES_INSTALL )
        install(
          FILES "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_RUNTIME}"
          DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
        install(
          FILES "$<TARGET_LINKER_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_IMPORT}"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
      endif()
    else()
      add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:${TARGET_NAME}>"
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_FILENAME}"
        VERBATIM
      )

      if ( SPLINES_INSTALL )
        install(
          FILES "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_FILENAME}"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
      endif()
    endif()
    return()
  endif()

  # Remove stale outputs before linking/archiving.  This prevents the old layout
  # (platform-suffixed symlink -> generic real library) from producing a
  # symlink cycle when switching to the new layout.
  add_custom_command(
    TARGET ${TARGET_NAME} PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E rm -f
            "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_FILENAME}"
    COMMAND ${CMAKE_COMMAND} -E rm -f
            "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_REAL_FILENAME}"
    VERBATIM
  )

  add_custom_command(
    TARGET ${TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f
            "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_FILENAME}"
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${TARGET_REAL_FILENAME}"
            "$<TARGET_FILE_DIR:${TARGET_NAME}>/${TARGET_LINK_FILENAME}"
    VERBATIM
  )

  if ( SPLINES_INSTALL )
    install(
      CODE
      "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
      execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f
              \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${TARGET_LINK_FILENAME}\")
      execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
              \"${TARGET_REAL_FILENAME}\"
              \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${TARGET_LINK_FILENAME}\")"
    )
  endif()
endfunction()

# Create the generic compatibility alias for the main Splines library in
# lib/lib after the real platform-suffixed artifact has been generated.
splines_add_platform_library_alias( Splines Splines )

#   _____         _
#  |_   _|__  ___| |_ ___
#    | |/ _ \/ __| __/ __|
#    | |  __/\__ \ |_\__ \
#    |_|\___||___/\__|___/
#
# Test executables.
#
# Tests are top-level only so a parent project that imports Splines with
# FetchContent does not inherit Splines' test battery.  Each test*.cc file under
# src_tests is converted into an executable and registered as a CTest test.
include( CTest )
if ( BUILD_TESTING AND PROJECT_IS_TOP_LEVEL )
  file(
    GLOB SPLINES_TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/test*.cc"
  )
  list( SORT SPLINES_TEST_SOURCES )
  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )
  foreach( TEST_SOURCE ${SPLINES_TEST_SOURCES} )
    get_filename_component( S "${TEST_SOURCE}" NAME_WE )
    add_executable( ${S} "${TEST_SOURCE}" )
    target_link_libraries( ${S} PRIVATE Splines::Splines ${CMAKE_DL_LIBS} )
    set_target_properties(
      ${S} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    add_test( NAME "${S}" COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${S}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
  endforeach()
endif()

#   ___                _                         _
#  | _ ) ___ _ _  __ | |_  _ __  __ _ _ _ | |__ ___
#  | _ \/ -_) ' \/ _|| ' \| '  \/ _` | '_|| / /(_-<
#  |___/\___|_||_\__||_||_|_|_|_\__,_|_|  |_\_\/__/
#
# Optional benchmark target.
#
# Benchmarks are intentionally top-level only, matching the test behavior.  This
# avoids polluting dependency builds with performance binaries.
if ( SPLINES_BUILD_BENCHMARKS AND PROJECT_IS_TOP_LEVEL )
  add_executable( bench_eval "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/bench_eval.cc" )
  target_link_libraries( bench_eval PRIVATE Splines::Splines )
  set_target_properties(
    bench_eval PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
# Install and package generation.
#
# The install target exports the Splines library and public headers into the
# configured local layout.  CPack is enabled only when installation is enabled,
# so dependency builds do not automatically gain package-generation targets.
if ( SPLINES_INSTALL )
  install(
    TARGETS Splines
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  set( CPACK_PACKAGE_NAME "${PROJECT_NAME}" )
  set( CPACK_PACKAGE_VERSION "${PROJECT_VERSION}" )
  set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}" )
  set( CPACK_GENERATOR "TGZ;ZIP" )
  set( CPACK_SOURCE_GENERATOR "TGZ;ZIP" )
  set( CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/" )
  include( CPack )
endif()

# ============================================================================
# Configuration summary
# ============================================================================
# Print the resolved dependency sources at the end of configuration.  These
# messages are useful for CI logs and for diagnosing whether a build used local
# sibling checkouts or downloaded dependencies through FetchContent.
message( STATUS "" )
message( STATUS "Splines dependency summary:" )
message( STATUS "  Eigen             : ${SPLINES_EIGEN_SOURCE}" )
message( STATUS "  GenericContainer  : ${GenericContainer_SOURCE}" )
message( STATUS "  UtilsLite         : ${UtilsLite_SOURCE}" )
message( STATUS "  quarticRootsFlocke: ${quarticRootsFlocke_SOURCE}" )
message( STATUS "  nlohmann_json     : ${SPLINES_NLOHMANN_JSON_SOURCE}" )
message( STATUS "" )
