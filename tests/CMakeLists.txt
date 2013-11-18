
# tests with CPU backend
foreach(PROG poisson_2d poisson_3d
            )
  add_executable(${PROG}-test src/${PROG}.cpp)
  add_test(${PROG} ${PROG}-test)
endforeach(PROG)

ADD_EXECUTABLE(external_linkage-test  src/external_1.cpp src/external_2.cpp)
ADD_test(external_linkage external_linkage-test)

