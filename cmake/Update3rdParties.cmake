############################################################################
#                                                                          #
#  file: cmake/Update3rdParties.cmake                                      #
#                                                                          #
#  Refresh third-party headers and stage MATLAB toolbox sources.          #
#                                                                          #
#  This replaces the old collect_dependencies post-build target with a     #
#  deterministic configure-time sync based on the dependency roots that    #
#  the top-level CMakeLists.txt has already resolved (local sibling        #
#  checkout first, FetchContent fallback otherwise).                       #
#                                                                          #
############################################################################

include_guard(GLOBAL)
include(CMakeParseArguments)

function(_splines_reset_dir DST)
  file(REMOVE_RECURSE "${DST}")
  file(MAKE_DIRECTORY "${DST}")
endfunction()

function(_splines_copy_matching_files SRC DST)
  if(NOT EXISTS "${SRC}")
    return()
  endif()

  foreach(_pattern IN LISTS ARGN)
    file(GLOB_RECURSE _matches RELATIVE "${SRC}" "${SRC}/${_pattern}")
    foreach(_rel IN LISTS _matches)
      if(_rel MATCHES "(^|/)(\\.DS_Store|\\.git|\\.github|\\.vscode)(/|$)")
        continue()
      endif()
      if(_rel MATCHES "(^|/)(build|CMakeFiles)(/|$)")
        continue()
      endif()
      get_filename_component(_dir "${_rel}" DIRECTORY)
      file(MAKE_DIRECTORY "${DST}/${_dir}")
      file(COPY "${SRC}/${_rel}" DESTINATION "${DST}/${_dir}")
    endforeach()
  endforeach()
endfunction()

function(_splines_copy_header_tree SRC DST)
  _splines_copy_matching_files(
    "${SRC}" "${DST}"
    "*.h" "*.hh" "*.hpp" "*.hxx" "*.cxx"
  )
endfunction()

function(_splines_copy_public_source_tree SRC DST)
  _splines_copy_matching_files(
    "${SRC}" "${DST}"
    "*.h" "*.hh" "*.hpp" "*.hxx" "*.cxx" "*.cc"
  )
endfunction()

function(_splines_copy_utilslite_tree ROOT DST)
  if(EXISTS "${ROOT}/lib/include")
    message(STATUS "Vendoring UtilsLite headers from ${ROOT}/lib/include")
    _splines_copy_header_tree("${ROOT}/lib/include" "${DST}")
  elseif(EXISTS "${ROOT}/src")
    message(STATUS "Vendoring UtilsLite headers from ${ROOT}/src")
    _splines_copy_header_tree("${ROOT}/src" "${DST}")
  else()
    message(FATAL_ERROR "UtilsLite headers not found under ${ROOT}")
  endif()
endfunction()

function(_splines_copy_eigen_tree ROOT DST)
  if(NOT EXISTS "${ROOT}/Eigen")
    message(FATAL_ERROR "Eigen headers not found under ${ROOT}")
  endif()

  file(COPY "${ROOT}/Eigen" DESTINATION "${DST}")
  if(EXISTS "${ROOT}/unsupported")
    file(COPY "${ROOT}/unsupported" DESTINATION "${DST}")
  endif()
endfunction()

function(_splines_require_directory PATH LABEL)
  if(NOT IS_DIRECTORY "${PATH}")
    message(FATAL_ERROR "${LABEL} not found: ${PATH}")
  endif()
endfunction()

function(_splines_copy_complete_tree SRC DST LABEL)
  _splines_require_directory("${SRC}" "${LABEL}")
  file(MAKE_DIRECTORY "${DST}")
  file(COPY "${SRC}/" DESTINATION "${DST}")
endfunction()

function(_splines_replace_literal PATH OLD_VALUE NEW_VALUE)
  if(NOT EXISTS "${PATH}")
    message(FATAL_ERROR "Cannot patch missing toolbox source: ${PATH}")
  endif()

  file(READ "${PATH}" _contents)
  string(REPLACE "${OLD_VALUE}" "${NEW_VALUE}" _updated "${_contents}")
  if(NOT _updated STREQUAL _contents)
    file(WRITE "${PATH}" "${_updated}")
  endif()
endfunction()

# Build the source bundle consumed by toolbox/CMakeLists.txt.  The operation
# intentionally mirrors the former Ruby population behavior, but uses
# dependency roots already resolved by the top-level CMake build.  It therefore
# works with both local sibling checkouts and FetchContent downloads.
function(splines_populate_toolbox)
  set(_options)
  set(_one_value_args
    DESTINATION
    SPLINES_ROOT
    JSON_INCLUDE_DIR
    GENERIC_CONTAINER_ROOT
    UTILSLITE_ROOT
    QUARTIC_ROOTS_ROOT
  )
  cmake_parse_arguments(SPLINES_TOOLBOX "${_options}" "${_one_value_args}" "" ${ARGN})

  if(NOT SPLINES_TOOLBOX_DESTINATION)
    message(FATAL_ERROR "splines_populate_toolbox requires DESTINATION")
  endif()

  set(_dst "${SPLINES_TOOLBOX_DESTINATION}")
  get_filename_component(_toolbox_root "${_dst}" DIRECTORY)
  set(_src_mex_dir "${_toolbox_root}/src_mex")

  message(STATUS "==============================================================")
  message(STATUS "Populating MATLAB toolbox sources in ${_dst}")
  message(STATUS "==============================================================")

  _splines_reset_dir("${_dst}")
  file(MAKE_DIRECTORY "${_src_mex_dir}" "${_toolbox_root}/bin")

  # Remove stale MEX binaries just as the former Ruby population step did.
  file(GLOB _mex_outputs "${_toolbox_root}/bin/*.mex*")
  if(_mex_outputs)
    file(REMOVE ${_mex_outputs})
  endif()

  # Copy in the same order as the old script.  Later trees may intentionally
  # replace files with the same relative path from an earlier tree.
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_SPLINES_ROOT}/src" "${_dst}" "Splines sources")
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_SPLINES_ROOT}/include" "${_dst}" "Splines headers")
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_QUARTIC_ROOTS_ROOT}/src" "${_dst}" "quarticRootsFlocke sources")
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_UTILSLITE_ROOT}/src" "${_dst}" "UtilsLite sources")
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_GENERIC_CONTAINER_ROOT}/src" "${_dst}" "GenericContainer sources")
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_GENERIC_CONTAINER_ROOT}/include" "${_dst}" "GenericContainer headers")
  _splines_copy_complete_tree("${SPLINES_TOOLBOX_JSON_INCLUDE_DIR}" "${_dst}" "nlohmann_json headers")

  set(
    _gc_matlab_interface
    "${SPLINES_TOOLBOX_GENERIC_CONTAINER_ROOT}/matlab/GenericContainerInterface_matlab.cc"
  )
  if(NOT EXISTS "${_gc_matlab_interface}")
    message(FATAL_ERROR "GenericContainer MATLAB interface not found: ${_gc_matlab_interface}")
  endif()
  file(COPY "${_gc_matlab_interface}" DESTINATION "${_dst}")
  file(COPY "${_gc_matlab_interface}" DESTINATION "${_src_mex_dir}")

  foreach(_interface_file
    "${_dst}/GenericContainerInterface_matlab.cc"
    "${_src_mex_dir}/GenericContainerInterface_matlab.cc"
  )
    _splines_replace_literal("${_interface_file}" "GC_ASSERT(" "GC_assert(")
  endforeach()

  # MATLAB's bundle does not use the top-level Eigen copies, and these legacy
  # UtilsLite translation units must not be compiled into the toolbox library.
  file(REMOVE_RECURSE "${_dst}/Eigen" "${_dst}/unsupported")
  foreach(_legacy_source
    Utils_Poly.cc
    Utils_GG2D.cc
    Utils_HJPatternSearch.cc
    Utils_NelderMead.cc
    Utils_nonlinear_system_tests.cc
  )
    file(REMOVE "${_dst}/${_legacy_source}")
  endforeach()

  if(NOT EXISTS "${SPLINES_TOOLBOX_SPLINES_ROOT}/license.txt")
    message(FATAL_ERROR "Splines license not found: ${SPLINES_TOOLBOX_SPLINES_ROOT}/license.txt")
  endif()
  file(COPY "${SPLINES_TOOLBOX_SPLINES_ROOT}/license.txt" DESTINATION "${_toolbox_root}")

  message(STATUS "MATLAB toolbox sources populated without Ruby")
  message(STATUS "==============================================================")
endfunction()

function(splines_update_3rdparties)
  set(_options)
  set(_one_value_args
    DESTINATION
    EIGEN_ROOT
    JSON_INCLUDE_DIR
    GENERIC_CONTAINER_ROOT
    UTILSLITE_ROOT
    QUARTIC_ROOTS_ROOT
  )
  cmake_parse_arguments(SPLINES_3RD "${_options}" "${_one_value_args}" "" ${ARGN})

  if(NOT SPLINES_3RD_DESTINATION)
    message(FATAL_ERROR "splines_update_3rdparties requires DESTINATION")
  endif()

  set(_dst "${SPLINES_3RD_DESTINATION}")
  message(STATUS "==============================================================")
  message(STATUS "Refreshing Splines third-party headers in ${_dst}")
  message(STATUS "==============================================================")

  _splines_reset_dir("${_dst}")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib")

  message(STATUS "Vendoring Eigen headers from ${SPLINES_3RD_EIGEN_ROOT}")
  _splines_copy_eigen_tree("${SPLINES_3RD_EIGEN_ROOT}" "${_dst}")

  if(EXISTS "${SPLINES_3RD_JSON_INCLUDE_DIR}/nlohmann/json.hpp")
    message(STATUS "Vendoring nlohmann_json headers from ${SPLINES_3RD_JSON_INCLUDE_DIR}")
    _splines_copy_header_tree("${SPLINES_3RD_JSON_INCLUDE_DIR}" "${_dst}")
  else()
    message(FATAL_ERROR "nlohmann_json headers not found under ${SPLINES_3RD_JSON_INCLUDE_DIR}")
  endif()

  if(EXISTS "${SPLINES_3RD_GENERIC_CONTAINER_ROOT}/include")
    message(STATUS "Vendoring GenericContainer headers from ${SPLINES_3RD_GENERIC_CONTAINER_ROOT}/include")
    _splines_copy_header_tree("${SPLINES_3RD_GENERIC_CONTAINER_ROOT}/include" "${_dst}")
  else()
    message(FATAL_ERROR "GenericContainer headers not found under ${SPLINES_3RD_GENERIC_CONTAINER_ROOT}")
  endif()

  _splines_copy_utilslite_tree("${SPLINES_3RD_UTILSLITE_ROOT}" "${_dst}")

  if(EXISTS "${SPLINES_3RD_QUARTIC_ROOTS_ROOT}/src")
    message(STATUS "Vendoring quarticRootsFlocke public sources from ${SPLINES_3RD_QUARTIC_ROOTS_ROOT}/src")
    _splines_copy_public_source_tree("${SPLINES_3RD_QUARTIC_ROOTS_ROOT}/src" "${_dst}")
  else()
    message(FATAL_ERROR "quarticRootsFlocke sources not found under ${SPLINES_3RD_QUARTIC_ROOTS_ROOT}")
  endif()

  message(STATUS "==============================================================")
  message(STATUS "Splines third-party headers refreshed. Review with `git diff`.")
  message(STATUS "==============================================================")
endfunction()
