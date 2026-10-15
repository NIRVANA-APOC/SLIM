cmake_minimum_required(VERSION 3.15)
project(slim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops -march=native -flto -fno-fat-lto-objects -ffast-math")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
endif()

file(GLOB SLIM_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp)
list(FILTER SLIM_SOURCES EXCLUDE REGEX ".*_main\\.cpp$")

find_package(Threads REQUIRED)

add_library(slim STATIC ${SLIM_SOURCES})
target_include_directories(slim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(slim PUBLIC cxx_std_17)
target_link_libraries(slim PUBLIC Threads::Threads)

add_executable(slim_index src/slim_index_main.cpp)
target_link_libraries(slim_index PRIVATE slim)

add_executable(slim_query src/slim_query_main.cpp)
target_link_libraries(slim_query PRIVATE slim)

set(SLIM_BINARY_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(TARGET slim_index POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE_NAME:slim_index>
            ${CMAKE_CURRENT_BINARY_DIR}/slim_offline)
add_custom_command(TARGET slim_query POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE_NAME:slim_query>
            ${CMAKE_CURRENT_BINARY_DIR}/slim_online)
