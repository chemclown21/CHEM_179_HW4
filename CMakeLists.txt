cmake_minimum_required(VERSION 3.27)
project(CHEM_179_HW4)

set(CMAKE_CXX_STANDARD 23)

add_executable(CHEM_179_HW4 main.cpp)

include_directories(/Users/vittor/Downloads/armadillo-12.6.7/include)

#Find LAPACK
find_package(LAPACK REQUIRED)
find_library(ARMA_LIBRARY armadillo)
target_link_libraries(CHEM_179_HW4 ${ARMA_LIBRARY} ${LAPACK_LIBRARIES})

