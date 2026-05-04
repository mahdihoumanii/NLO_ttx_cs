cmake_minimum_required(VERSION 3.20)

project(cs_ppttb LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(LHAPDF_CONFIG_EXECUTABLE
    NAMES lhapdf-config
    HINTS
        ENV CONDA_PREFIX
        "$ENV{CONDA_PREFIX}/bin")

if(NOT LHAPDF_CONFIG_EXECUTABLE)
    message(FATAL_ERROR "Could not find lhapdf-config. Set LHAPDF_CONFIG_EXECUTABLE or activate the LHAPDF environment.")
endif()

execute_process(
    COMMAND ${LHAPDF_CONFIG_EXECUTABLE} --incdir
    OUTPUT_VARIABLE LHAPDF_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${LHAPDF_CONFIG_EXECUTABLE} --libdir
    OUTPUT_VARIABLE LHAPDF_LIBRARY_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

find_library(LHAPDF_LIBRARY
    NAMES LHAPDF
    PATHS ${LHAPDF_LIBRARY_DIR}
    NO_DEFAULT_PATH)

if(NOT LHAPDF_LIBRARY)
    message(FATAL_ERROR "Could not locate the LHAPDF library in ${LHAPDF_LIBRARY_DIR}")
endif()

add_executable(cs_ppttb
    apps/cs_ppttb.cpp
    src/born.cpp
    src/color_correlation.cpp
    src/virtual.cpp
    src/cs_born_correlations.cpp
    src/cs_ca_integrand.cpp
    src/cs_ct_shift.cpp
    src/cs_dipole_metadata.cpp
    src/cs_ioperator.cpp
    src/cs_kinematics.cpp
    src/cs_va_integrand.cpp
    src/cs_va_kernel.cpp
    src/cs_config.cpp
    src/cs_phasespace.cpp
    src/cs_plot.cpp
    src/cs_runner.cpp)

target_include_directories(cs_ppttb PRIVATE include ${LHAPDF_INCLUDE_DIR})
target_link_libraries(cs_ppttb PRIVATE ${LHAPDF_LIBRARY})

set_target_properties(cs_ppttb PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)