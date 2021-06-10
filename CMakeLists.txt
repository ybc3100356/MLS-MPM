cmake_minimum_required(VERSION 3.19)
project(mls-mpm)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenCV REQUIRED)
# directory of opencv headers
include_directories(${OpenCV_INCLUDE_DIRS})
# name of executable file and path of source file
add_executable(mls-mpm main.cpp)

# directory of opencv library
link_directories(${OpenCV_LIBRARY_DIRS})
# opencv libraries
target_link_libraries(mls-mpm ${OpenCV_LIBS})
