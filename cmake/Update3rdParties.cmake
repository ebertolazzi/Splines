############################################################################
#                                                                          #
#  file: cmake/Update3rdParties.cmake                                      #
#                                                                          #
#  Refresh the vendored third-party headers under lib3rd/include.          #
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

function(splines_update_3rdparties)
  set(_options)
  set(_one_value_args
    DESTINATION
    EIGEN_ROOT
    JSON_ROOT
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

  if(EXISTS "${SPLINES_3RD_JSON_ROOT}/include")
    message(STATUS "Vendoring nlohmann_json headers from ${SPLINES_3RD_JSON_ROOT}/include")
    _splines_copy_header_tree("${SPLINES_3RD_JSON_ROOT}/include" "${_dst}")
  else()
    message(FATAL_ERROR "nlohmann_json headers not found under ${SPLINES_3RD_JSON_ROOT}")
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
