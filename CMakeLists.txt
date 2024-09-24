cmake_minimum_required(VERSION 3.8)
project(sutoku)
include(ExternalProject)
option(INSTALL_SLIB "install slib" OFF) 
option(BUILD_SLIB "only build slib" OFF) 
if(MSVC)
	add_compile_definitions(_USE_DLL)
	set(CMAKE_CXX_FLAGS "/utf-8")
	link_libraries(sobjd.lib sbioinfod.lib sappd.lib)
elseif(UNIX)
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
	find_library(LSOBJ_INSTALLED NAMES sobj)
	find_library(LSBI_INSTALLED NAMES sbioinfo)
	find_library(LSAPP_INSTALLED NAMES sapp)
	if (LSOBJ_INSTALLED AND LSBI_INSTALLED AND LSAPP_INSTALLED)
		message(STATUS "slib found.")
	elseif (INSTALL_SLIB)
		message(STATUS "Error: Three libraries <libsobj, libsbioinfo, libsapp> are required.")
		ExternalProject_Add(
			slib
			PREFIX ${CMAKE_BINARY_DIR}
			GIT_REPOSITORY https://github.com/YujiSue/slib.git
			GIT_TAG master
    			CMAKE_ARGS -DBUILD_SHARED_LIBS=ON -DBUILD_BIOINFO_LIB=ON -DBUILD_SAPP_LIB=ON 
    			UPDATE_DISCONNECTED 1
		)
	elseif (BUILD_SLIB)
		message(STATUS "Error: Three libraries <libsobj, libsbioinfo, libsapp> are required.")
		ExternalProject_Add(
			slib
			PREFIX ${CMAKE_SOURCE_DIR}
			GIT_REPOSITORY https://github.com/YujiSue/slib.git
			GIT_TAG master
    			CMAKE_ARGS -DBUILD_ONLY=ON -DBUILD_SHARED_LIBS=ON -DBUILD_BIOINFO_LIB=ON -DBUILD_SAPP_LIB=ON 
    			UPDATE_DISCONNECTED 1
		)
		ExternalProject_Get_Property(slib install_dir)
		message(STATUS "slib install : ${install_dir}")
		include_directories(${install_dir}/src/slib/include)
		link_directories(${install_dir}/src/slib-build/lib)	
	endif()
	message(STATUS "Install to ${CMAKE_INSTALL_PREFIX}")
	set(CMAKE_CXX_FLAGS "-std=c++11 -O2")
	link_libraries(sobj sbioinfo sapp dl pthread curl)
endif()
include_directories(${PROJECT_SOURCE_DIR} $ENV{SLIB_INCLUDE})
link_directories($ENV{SLIB_DIR})
execute_process(
	COMMAND python MakeProf.py profile.h.prv appinfo.json ./src/profile.h
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_subdirectory(src)
