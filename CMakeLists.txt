# set up project
project("dune-pdelab-howto" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

# guess dune-common build dir
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

# we need the module file to be able to build via dunecontrol
set(EXTRA_DIST             
        dune.module     
        BSD             
        CC-BY-SA        
        RELEASE_NOTES)

set(SUBDIRS  doc m4 src)

# don't follow the full GNU-standard
# we need automake 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)

# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/top-rules
# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})

add_subdirectory(cmake/modules)
# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
