cmake_minimum_required(VERSION 3.2)
project(c4cfd)

include(./cmake/c4Project.cmake)
c4_declare_project(c4cfd)

set(C4CFD_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
set(C4CFD_EXT_DIR ${CMAKE_CURRENT_LIST_DIR}/ext)

c4_require_subproject(c4core REMOTE
    GIT_REPOSITORY https://github.com/biojppm/c4core.git)
c4_require_subproject(ryml REMOTE
    GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git)

c4_add_library(c4cfd
    SOURCE_ROOT ${C4CFD_SRC_DIR}
    SOURCES
        c4/cfd/cfd.hpp
        c4/cfd/cfd.cpp
    LIBS c4core
    INC_DIRS
       $<BUILD_INTERFACE:${C4CFD_SRC_DIR}> $<INSTALL_INTERFACE:include>
)

c4_install_target(c4cfd)
c4_install_exports()

c4_add_dev_targets()

