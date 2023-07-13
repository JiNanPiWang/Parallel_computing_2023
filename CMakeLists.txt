cmake_minimum_required(VERSION 3.25)
project(并行计算)

# 添加子目录
add_subdirectory(codes_1)
add_subdirectory(codes_2)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
