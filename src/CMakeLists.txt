cmake_minimum_required(VERSION 3.10)
project(HTS2D)

target_include_directories(HTS2D /work/ese-chenph/lib/petsc-3.16.6-opt/include)
add_executable(HTS2d main.c explicit.c implicit.c plicit.h)
target_link_libraries(HTS2d /work/ese-chenph/lib/petsc-3.16.6-opt/lib)