##
# PYTHON Wrapper for Splines
#
# License MIT - See LICENSE file
# 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
#

cmake_minimum_required(VERSION 3.28)

project(Splines VERSION 1.0 LANGUAGES CXX)

# -----------------------------------------------------------------------------
# C++ standard
# Utils in lib3rd requires C++17.
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# Python selection
#
# With modern pybind11, PYTHON_EXECUTABLE is ignored because pybind11 uses
# CMake's FindPython backend. Use Python_EXECUTABLE instead.
#
# Recommended configure command:
#   cmake -S . -B build -DPython_EXECUTABLE:FILEPATH="$CONDA_PREFIX/bin/python"
# -----------------------------------------------------------------------------
set(PYBIND11_FINDPYTHON ON CACHE BOOL "Use CMake FindPython backend for pybind11")

if(DEFINED PYTHON_EXECUTABLE AND NOT DEFINED Python_EXECUTABLE)
  set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}" CACHE FILEPATH "Python executable" FORCE)
endif()

if(DEFINED ENV{CONDA_PREFIX} AND NOT DEFINED Python_EXECUTABLE)
  set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python" CACHE FILEPATH "Python executable" FORCE)
endif()

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
message(STATUS "Python executable: ${Python_EXECUTABLE}")
message(STATUS "Python version:    ${Python_VERSION}")

find_package(pybind11 REQUIRED)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
set(SPLINES "${CMAKE_CURRENT_SOURCE_DIR}/../lib")
set(GENERIC_CONTAINER "${CMAKE_CURRENT_SOURCE_DIR}/../lib3rd")

set(SOURCES
  "src/py-Splines.cc"
  "src/py-CubicSplines.cc"
  "src/py-ConstantSplines.cc"
  "src/py-LinearSplines.cc"
  "src/py-AkimaSplines.cc"
  "src/py-VanLeerSplines.cc"
  "src/py-HermiteSplines.cc"
  "src/py-PchipSplines.cc"
  "src/py-QuinticSplines.cc")

pybind11_add_module(${PROJECT_NAME} MODULE ${SOURCES} NO_EXTRAS)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON)

target_include_directories(${PROJECT_NAME}
  PUBLIC "${SPLINES}/include"
  PUBLIC "${GENERIC_CONTAINER}/include"
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(APPLE)
  set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_VISIBILITY_PRESET default)
  set(GENERIC_CONTAINER_LIB "${GENERIC_CONTAINER}/lib/libGenericContainer_osx_static.a")
  set(SPLINES_LIB "${SPLINES}/lib/libSplines_osx_static.a")
  set(QUARTIC_ROOTS_LIB "${GENERIC_CONTAINER}/lib/libquarticRootsFlocke_osx_static.a")
elseif(UNIX)
  set(GENERIC_CONTAINER_LIB "${GENERIC_CONTAINER}/lib/libGenericContainer_linux_static.a")
  set(SPLINES_LIB "${SPLINES}/lib/libSplines_linux_static.a")
  set(QUARTIC_ROOTS_LIB "${GENERIC_CONTAINER}/lib/libquarticRootsFlocke_linux_static.a")
else()
  set(GENERIC_CONTAINER_LIB "${GENERIC_CONTAINER}/lib/GenericContainer_win__static.lib")
  set(SPLINES_LIB "${SPLINES}/lib/Splines_win__static.lib")
  set(QUARTIC_ROOTS_LIB "${GENERIC_CONTAINER}/lib/quarticRootsFlocke_win__static.lib")
endif()

if(NOT EXISTS "${GENERIC_CONTAINER_LIB}")
  message(FATAL_ERROR
    "GenericContainer library not found: ${GENERIC_CONTAINER_LIB}\n"
    "Expected headers in: ${GENERIC_CONTAINER}/include\n"
    "Expected library in: ${GENERIC_CONTAINER}/lib")
endif()

if(NOT EXISTS "${SPLINES_LIB}")
  message(FATAL_ERROR "Splines library not found: ${SPLINES_LIB}")
endif()

if(NOT EXISTS "${QUARTIC_ROOTS_LIB}")
  message(FATAL_ERROR "quarticRootsFlocke library not found: ${QUARTIC_ROOTS_LIB}")
endif()

target_link_libraries(${PROJECT_NAME}
  PUBLIC "${GENERIC_CONTAINER_LIB}"
  PUBLIC "${SPLINES_LIB}"
  PUBLIC "${QUARTIC_ROOTS_LIB}")

# -----------------------------------------------------------------------------
# Install
# -----------------------------------------------------------------------------
set(INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../lib/lib")
install(TARGETS ${PROJECT_NAME} DESTINATION "${INSTALL_PREFIX}")

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
enable_testing()

set(TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/test_splines.py")
if(EXISTS "${TEST_SCRIPT}")
  add_test(
    NAME python_import_splines
    COMMAND "${Python_EXECUTABLE}" "${TEST_SCRIPT}" --module-dir "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
  )

  set_tests_properties(python_import_splines PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PROJECT_NAME}>;${INSTALL_PREFIX};$ENV{PYTHONPATH}"
  )

  add_custom_target(tests
    COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  )
endif()
