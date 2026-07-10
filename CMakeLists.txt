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

# Estrai la versione da Git tag (se disponibile)
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
# Determinazione del suffisso per la piattaforma (come in GenericContainer)
# ============================================================================
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
# Opzioni e impostazioni di base
# ============================================================================

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
set( SPLINES_GC_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../GC" CACHE PATH "Local GenericContainer checkout" )
set( SPLINES_UTILSLITE_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite" CACHE PATH "Local UtilsLite checkout" )
set( SPLINES_QUARTIC_ROOTS_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../quarticRootsFlocke" CACHE PATH "Local quarticRootsFlocke checkout" )
set( SPLINES_EIGEN_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Eigen" CACHE PATH "Local Eigen checkout" )
set( SPLINES_JSON_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../json" CACHE PATH "Local nlohmann_json checkout" )
set( SPLINES_EIGEN_VERSION "5.0.0" CACHE STRING "Eigen version to fetch" )

# Default install layout requested for this project:
#   ${CMAKE_CURRENT_SOURCE_DIR}/lib/include
#   ${CMAKE_CURRENT_SOURCE_DIR}/lib/lib
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE )
endif()
set( CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE )
set( CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE )

# NOTA: non impostiamo più CMAKE_*_OUTPUT_DIRECTORY globalmente;
#       le directory di output verranno assegnate per singolo target.

include( GNUInstallDirs )
include( FetchContent )
include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Update3rdParties.cmake" )

if ( APPLE )
  execute_process( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
  set( CMAKE_OSX_ARCHITECTURES ${ARCHITECTURE} )
endif()

# ============================================================================
# Dependency helpers
# ============================================================================
function( splines_use_local_or_fetch NAME LOCAL_DIR )
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  string( TOUPPER "${NAME}" _name_upper )
  string( TOLOWER "${NAME}" _name_lower )
  if ( FETCHCONTENT_BASE_DIR )
    set( _fetchcontent_base_dir "${FETCHCONTENT_BASE_DIR}" )
  else()
    set( _fetchcontent_base_dir "${CMAKE_BINARY_DIR}/_deps" )
  endif()
  set( _dependency_install_stub "${_fetchcontent_base_dir}/${_name_lower}-build/cmake_install.cmake" )
  set( _saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}" )
  set( CMAKE_SKIP_INSTALL_RULES ON )

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
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( nlohmann_json "${SPLINES_JSON_LOCAL_DIR}" )
set( SPLINES_NLOHMANN_JSON_SOURCE "${nlohmann_json_SOURCE}" )



# ============================================================================
# Eigen
# ============================================================================
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
# GenericContainer
# ============================================================================
set( GENERIC_CONTAINER_INSTALL OFF CACHE BOOL "Disable GenericContainer install rules inside Splines" FORCE )
FetchContent_Declare(
  GenericContainer
  GIT_REPOSITORY https://github.com/ebertolazzi/GenericContainer.git
  GIT_TAG        GC2.0
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( GenericContainer "${SPLINES_GC_LOCAL_DIR}" )
set( GC GenericContainer::Yaml )

# Patch per Clang: forzare -stdlib=libc++
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

# UtilsLite and quarticRootsFlocke share the same vendored cmake_utils option
# to gate their own internal test suites -- turn it off so ctest here only
# runs Splines' own test00..test18.
set( UTILS_INSTALL_LOCAL OFF CACHE BOOL "Disable dependency-local install rules inside Splines" FORCE )
set( UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable dependency test suites" FORCE )

# ============================================================================
# UtilsLite
# ============================================================================
FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        Utils2.0
  GIT_SHALLOW    TRUE
)
splines_use_local_or_fetch( UtilsLite "${SPLINES_UTILSLITE_LOCAL_DIR}" )
set( UTILS utils::UtilsLite_Static )

# Treat UtilsLite's headers as SYSTEM includes for anything that links it.
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
# quarticRootsFlocke
# ============================================================================
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
# Funzione per impostare le directory di output per un target
# ============================================================================
function( splines_set_target_output_dirs TARGET_NAME LIB_DIR )
  set_target_properties( ${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

# ============================================================================
# Funzione per impostare le directory di output per le dipendenze
# (gestisce alias target)
# ============================================================================
function( splines_set_dependency_output_dirs LIB_DIR )
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

# Imposta le directory di output per Splines (vanno in lib/lib)
splines_set_target_output_dirs( Splines lib/lib )

target_compile_features( Splines PUBLIC cxx_std_20 )
target_compile_definitions( Splines PUBLIC SPLINES_NO_COMPATIBILITY )
target_compile_definitions( Splines PUBLIC _USE_MATH_DEFINES )
target_include_directories(
  Splines PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries( Splines PUBLIC ${ROOTS} ${GC} ${UTILS} nlohmann_json::nlohmann_json )

set_target_properties(
  Splines PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Splines
)

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

if ( MSVC )
  target_compile_options( Splines PUBLIC /utf-8 /bigobj )
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )

# ============================================================================
# Imposta le directory di output per le dipendenze (vanno in lib3rd/lib)
# ============================================================================
splines_set_dependency_output_dirs( lib3rd/lib
  GenericContainer
  Lua
  ${UTILS}
  ${ROOTS}
)

# ============================================================================
# Creazione dell'alias con suffisso OS (come in GenericContainer)
# ============================================================================
function( splines_add_platform_library_alias TARGET_NAME BASE_NAME )
  set( ALIAS_BASE "lib${BASE_NAME}_${SPLINES_PLATFORM_SUFFIX}" )
  set( LEGACY_ALIAS_BASE "${ALIAS_BASE}" )

  if ( NOT BUILD_SHARED_LIBS )
    string( APPEND ALIAS_BASE "_static" )
  endif()

  if ( BUILD_SHARED_LIBS )
    set( LINK_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
    set( ALIAS_FILENAME "${ALIAS_BASE}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
  else()
    set( LINK_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    set( ALIAS_FILENAME "${ALIAS_BASE}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    set( LEGACY_ALIAS_FILENAME "${LEGACY_ALIAS_BASE}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
  endif()

  if ( WIN32 )
    if ( BUILD_SHARED_LIBS )
      add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:${TARGET_NAME}>"
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
        VERBATIM
      )
    else()
      add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E rm -f
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${LEGACY_ALIAS_FILENAME}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:${TARGET_NAME}>"
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
        VERBATIM
      )
      install(
        CODE
        "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${LEGACY_ALIAS_FILENAME}\")"
      )
    endif()
    install(
      FILES "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  else()
    if ( BUILD_SHARED_LIBS )
      add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E rm -f "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${LINK_FILENAME}"
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
        VERBATIM
      )
    else()
      add_custom_command(
        TARGET ${TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E rm -f "$<TARGET_FILE_DIR:${TARGET_NAME}>/${LEGACY_ALIAS_FILENAME}"
        COMMAND ${CMAKE_COMMAND} -E rm -f "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${LINK_FILENAME}"
                "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ALIAS_FILENAME}"
        VERBATIM
      )
    endif()
    if ( BUILD_SHARED_LIBS )
      install(
        CODE
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
        execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\")
        execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
                \"${LINK_FILENAME}\"
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\")"
      )
    else()
      install(
        CODE
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
        execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${LEGACY_ALIAS_FILENAME}\")
        execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\")
        execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
                \"${LINK_FILENAME}\"
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${ALIAS_FILENAME}\")"
      )
    endif()
  endif()
endfunction()

# Crea alias per Splines (in lib/lib)
splines_add_platform_library_alias( Splines Splines )

#   _____         _
#  |_   _|__  ___| |_ ___
#    | |/ _ \/ __| __/ __|
#    | |  __/\__ \ |_\__ \
#    |_|\___||___/\__|___/
#
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
# Riepilogo delle dipendenze
# ============================================================================
message( STATUS "" )
message( STATUS "Splines dependency summary:" )
message( STATUS "  Eigen             : ${SPLINES_EIGEN_SOURCE}" )
message( STATUS "  GenericContainer  : ${GenericContainer_SOURCE}" )
message( STATUS "  UtilsLite         : ${UtilsLite_SOURCE}" )
message( STATUS "  quarticRootsFlocke: ${quarticRootsFlocke_SOURCE}" )
message( STATUS "  nlohmann_json     : ${SPLINES_NLOHMANN_JSON_SOURCE}" )
message( STATUS "" )
