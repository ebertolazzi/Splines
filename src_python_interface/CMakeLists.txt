##
# PYTHON Wrapper for Clothoids
# 
# License MIT - See LICENSE file
# 2019 Matteo Ragni, Claudio Kerov Ghiglianovich, Enrico Bertolazzi
#

cmake_minimum_required(VERSION 3.10.2)

project(Splines VERSION 1.0)

find_package(pybind11 REQUIRED)

find_package(pybind11)

	set(SPLINES "${CMAKE_CURRENT_SOURCE_DIR}/../lib")
	set(GENERIC_CONTAINER "${CMAKE_CURRENT_SOURCE_DIR}/../GC/lib")

	set(SOURCES
		"src/py-Splines.cc"
		"src/py-CubicSplines.cc"
		"src/py-ConstantSplines.cc"
		"src/py-LinearSplines.cc"
		"src/py-AkimaSplines.cc"
		"src/py-BesselSplines.cc"
		"src/py-HermiteSplines.cc"
		"src/py-PchipSplines.cc"
		"src/py-QuinticSplines.cc")

pybind11_add_module(${PROJECT_NAME} MODULE ${SOURCES} NO_EXTRAS)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

target_include_directories(${PROJECT_NAME}
	PUBLIC "${SPLINES}/include"
	PUBLIC "${GENERIC_CONTAINER}/include"
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

set_property(TARGET ${PROJECT_NAME} PROPERTY
		POSITION_INDEPENDENT_CODE 1)

if(APPLE)
  set_property(TARGET ${PROJECT_NAME}
    PROPERTY CXX_VISIBILITY_PRESET default)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC "${GENERIC_CONTAINER}/lib/libGenericContainer_osx_static.a")
  target_link_libraries(${PROJECT_NAME}
    PUBLIC "${SPLINES}/lib/libSplines_osx_static.a")
ELSEIF(LINUX OR UNIX)
  target_link_libraries(${PROJECT_NAME}
    PUBLIC "${GENERIC_CONTAINER}/lib/libGenericContainer_linux_static.a")
  target_link_libraries(${PROJECT_NAME}
    PUBLIC "${SPLINES}/lib/libSplines_linux_static.a")
else()
  target_link_libraries(${PROJECT_NAME}
    PUBLIC "${GENERIC_CONTAINER}/lib/GenericContainer_win__static.lib")
  target_link_libraries(${PROJECT_NAME}
    PUBLIC "${SPLINES}/lib/Splines_win__static.lib")
endif()

set(INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/../lib/lib)
install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_PREFIX})
