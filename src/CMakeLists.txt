set(src
    c4/cfd/cfd.hpp
    c4/cfd/cfd.cpp
)

c4_require_module(c4cfd
    c4core REMOTE GIT_REPOSITORY https://github.com/biojppm/c4core.git)
#c4_require_module(c4cfd
#    rapidyaml REMOTE GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git)

c4_add_library(C4CFD c4cfd
    SOURCES ${src}
    LIBS c4core
    INC_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
)
