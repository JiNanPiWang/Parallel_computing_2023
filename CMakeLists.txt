cmake_minimum_required(VERSION 3.25)
set(CMAKE_C_STANDARD 11)
project(并行计算)

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

# 添加子目录
add_subdirectory(codes_1)
add_subdirectory(codes_2)
add_subdirectory(codes_3)
add_subdirectory(codes_4)
add_subdirectory(codes_5)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
