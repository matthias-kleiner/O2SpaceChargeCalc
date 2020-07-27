cmake_minimum_required(VERSION 3.12)
project(MyPowerProject CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ standard: C++${CMAKE_CXX_STANDARD}")

find_package(ROOT)
if(ROOT_FOUND)
    include(${ROOT_USE_FILE} )
    message(STATUS "Using ROOT: ${ROOT_VERSION} <${ROOT_CONFIG}>")
endif(ROOT_FOUND)

# INCLUDE_DIRECTORIES ( "/Users/matthias/Desktop/testCube/eigen-3.3.7/Eigen/" )

set(CMAKE_CXX_FLAGS_DEBUG "-I${VC_ROOT}/include")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -march=native")

set(CMAKE_CXX_FLAGS_DEBUG  "-g -O3")

add_executable(executable main.cpp include/TriCubic.h include/RegularGrid3D.h include/matrix.h include/O2TPCSpaceCharge3DCalc.h include/AliRoot/AliTPCPoissonSolver.cxx)
find_package(OpenMP REQUIRED)
find_package(Vc REQUIRED)


target_link_libraries(executable PRIVATE OpenMP::OpenMP_CXX Vc::Vc ${ROOT_LIBRARIES} Vc::Vc "${TORCH_LIBRARIES}")
