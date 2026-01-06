cmake_minimum_required(VERSION 3.15)
project(slim VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find dependencies
find_package(OpenMP REQUIRED)

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mavx2 -DNDEBUG -flto -fno-fat-lto-objects -funroll-loops -ffast-math")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto -O3 -s")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")

if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
endif()

# Optional allocators
find_library(TCMALLOC_LIB tcmalloc)
find_library(JEMALLOC_LIB jemalloc)

# Add subdirectories
add_subdirectory(common)
add_subdirectory(online)
add_subdirectory(offline)
