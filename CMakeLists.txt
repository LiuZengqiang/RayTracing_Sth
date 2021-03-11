project(RayTracing_Sth)
cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 14)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)


add_executable(RayTracing_Sth main.cpp src/Scene.cpp)