cmake_minimum_required(VERSION 3.6)
project(afsr)

set(CMAKE_CXX_STANDARD 11)

#Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})

#Boost
find_package( Boost 1.68.0 REQUIRED COMPONENTS system filesystem regex)
include_directories( ${Boost_INCLUDE_DIRS})

#CGAL
set(CGAL_DIR "/usr/local/")

find_package(CGAL COMPONENTS Core REQUIRED)
include(${CGAL_USE_FILE})
include( CGAL_CreateSingleSourceCGALProgram )
message("CGAL DIR ${CGAL_DIR}")

add_executable(afsr main.cpp)
target_link_libraries(afsr
        PRIVATE
        Boost::system
        Boost::filesystem
        Boost::regex
        )
