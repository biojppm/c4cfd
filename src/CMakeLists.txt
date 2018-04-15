set(src
    c4/cfd/cfd.hpp
    c4/cfd/cfd.cpp
)

c4_add_library(C4CFD c4cfd
    SOURCES ${src}
    LIBS c4core
    INC_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
)
