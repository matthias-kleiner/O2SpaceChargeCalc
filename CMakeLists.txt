cmake_minimum_required(VERSION 3.4...3.16)
project(SpaceCharge3D CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ standard: C++${CMAKE_CXX_STANDARD}")

find_package(ROOT)
if(ROOT_FOUND)
    include(${ROOT_USE_FILE} )
    message(STATUS "Using ROOT: ${ROOT_VERSION} <${ROOT_CONFIG}>")
endif(ROOT_FOUND)

#find_package(ROOT 6.23 CONFIG REQUIRED)


#set(CMAKE_CXX_FLAGS "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wall -Wextra -O3 -march=native -pedantic")

set(CMAKE_CXX_FLAGS "-Wall -O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG  "-g -O3")

find_package(OpenMP REQUIRED)
find_package(Vc REQUIRED)
find_package(CGAL REQUIRED)

link_libraries(Vc::Vc)
link_libraries(ROOT::Core)

ROOT_GENERATE_DICTIONARY(G__RegularGrid3D src/include/RegularGrid3D.h LINKDEF src/include/SpaceChargeLinkDef.h)
add_library(RegularGrid3D SHARED src/include/RegularGrid3D.cxx src/include/RegularGrid3D.h G__RegularGrid3D.cxx)
target_include_directories(RegularGrid3D PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(executable main.cpp src/include/TriCubic.h src/include/RegularGrid3D.h src/include/Vector.h src/O2TPCSpaceCharge3DCalc.cxx src/include/AliRoot/AliTPCPoissonSolver.cxx)

target_link_libraries(executable PRIVATE RegularGrid3D OpenMP::OpenMP_CXX Vc::Vc ROOT::Core ROOT::Physics ROOT::ROOTDataFrame)
