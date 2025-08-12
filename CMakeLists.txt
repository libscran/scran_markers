cmake_minimum_required(VERSION 3.14)

project(scran_markers
    VERSION 0.1.4
    DESCRIPTION "Marker detection for single-cell data"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(scran_markers INTERFACE)
add_library(libscran::scran_markers ALIAS scran_markers)

target_include_directories(scran_markers INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scran_markers>)
target_compile_features(scran_markers INTERFACE cxx_std_17)

# Dependencies
option(SCRAN_MARKERS_FETCH_EXTERN "Automatically fetch scran_markers's external dependencies." ON)
if(SCRAN_MARKERS_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 4.0.0 CONFIG REQUIRED)
    find_package(tatami_tatami_stats 2.0.0 CONFIG REQUIRED)
    find_package(libscran_scran_blocks 0.1.0 CONFIG REQUIRED)
    find_package(ltla_sanisizer 0.1.1 CONFIG REQUIRED)
endif()

target_link_libraries(scran_markers INTERFACE tatami::tatami tatami::tatami_stats libscran::scran_blocks ltla::sanisizer)

# Tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(SCRAN_MARKERS_TESTS "Build scran_markers's test suite." ON)
else()
    option(SCRAN_MARKERS_TESTS "Build scran_markers's test suite." OFF)
endif()

if(SCRAN_MARKERS_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/scran_markers)

install(TARGETS scran_markers
    EXPORT scran_markersTargets)

install(EXPORT scran_markersTargets
    FILE libscran_scran_markersTargets.cmake
    NAMESPACE libscran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_markers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_markersConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_markers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_markersConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_markersConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_markersConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_markers)
