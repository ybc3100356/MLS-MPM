cmake_minimum_required(VERSION 3.19)
project(mls-mpm)
enable_language(CUDA)

set(CMAKE_CXX_STANDARD 14)
find_package(glfw3 3.3 REQUIRED)

file(GLOB SOURCES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/src/*.c
        ${PROJECT_SOURCE_DIR}/src/*.cu
        ${PROJECT_SOURCE_DIR}/src/shader/*.cpp
        )

add_executable(mls-mpm ${SOURCES})
set_target_properties(mls-mpm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#target_compile_features(mls-mpm PUBLIC cxx_std_11)

target_include_directories(mls-mpm PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(mls-mpm glfw)
