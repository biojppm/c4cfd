c4_setup_testing()

function(c4cfd_add_test name)
    c4_add_executable(c4cfd-test-${name}
        SOURCES ${ARGN} main.cpp
        LIBS c4cfd gtest
        FOLDER test)
    c4_add_test(c4cfd-test-${name})
endfunction(c4cfd_add_test)

c4cfd_add_test(basic basic.cpp)
