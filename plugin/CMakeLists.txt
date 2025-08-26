cmake_minimum_required(VERSION 3.9)
project(sutokuplugin)

set(CMAKE_SHARED_LIBRARY_PREFIX "")
set(CMAKE_SHARED_LIBRARY_SUFFIX ".plugin")

if(MSVC)
    add_compile_definitions(_USE_DLL)
    set(CMAKE_CXX_FLAGS_INIT "/utf-8")
    link_libraries(sobj.lib sbioinfo.lib sapp.lib)
else()
    if (DEFINED ENV{SLIB_INCLUDE})
        message(STATUS "SLIB_INCLUDE=$ENV{SLIB_INCLUDE}")
    else()
        set(ENV{SLIB_INCLUDE} ${CMAKE_INSTALL_PREFIX}/include/slib)
        message(STATUS "SLIB_INCLUDE=$ENV{SLIB_INCLUDE}")
    endif()
    if (DEFINED ENV{SLIB_DIR})
        message(STATUS "SLIB_DIR=$ENV{SLIB_DIR}")
    else()
        set(ENV{SLIB_DIR} ${CMAKE_INSTALL_PREFIX}/lib)
        message(STATUS "SLIB_DIR=$ENV{SLIB_DIR}")
    endif()
    set(CMAKE_CXX_FLAGS_INIT "-std=c++11 -O2")
    link_libraries(sobj sbioinfo sapp)
endif()

include_directories("${PROJECT_SOURCE_DIR}/src" $ENV{SLIB_INCLUDE})
link_directories($ENV{SLIB_DIR})

# MutagenHomoFilter.cpp
add_library(MutagenHomoFilter SHARED MutagenHomoFilter.cpp)

# VCHardFilter.cpp
add_library(VCHardFilter SHARED VCHardFilter.cpp)
