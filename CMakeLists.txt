cmake_minimum_required(VERSION 3.40)
project(PathTracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(path-tracer
    src/main.cpp
)

target_link_libraries(path-tracer PRIVATE OpenMP::OpenMP_CXX)