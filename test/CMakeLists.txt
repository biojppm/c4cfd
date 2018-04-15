c4_setup_testing(C4CFD ON)

function(c4cfd_add_test name)
    c4_add_executable(C4CFD c4cfd-test-${name}
        SOURCES ${ARGN} main.cpp
        INC_DIRS ${CMAKE_CURRENT_LIST_DIR}
        LIBS c4cfd gtest
        FOLDER c4cfd-test)
    c4_add_test(C4CFD c4cfd-test-${name} ON)
endfunction(c4cfd_add_test)

c4cfd_add_test(basic     basic.cpp)
