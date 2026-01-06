set(ONLINE_SOURCES
    src/filter.cpp
    src/matcher.cpp
    src/enumerate/enumerator.cpp
    src/enumerate/candidate_storage.cpp
    src/enumerate/unified_edge_matrix.cpp
    src/enumerate/unified_enumerator.cpp
    src/intersection/avx2.cpp
    src/intersection/standard.cpp
)

set(ONLINE_HEADERS
    include/slim/filter.h
    include/slim/matcher.h
    include/slim/enumerate/enumerator.h
    include/slim/enumerate/candidate_storage.h
    include/slim/enumerate/unified_edge_matrix.h
    include/slim/intersection/avx2.h
    include/slim/intersection/standard.h
)

add_library(slim_online_lib STATIC ${ONLINE_SOURCES} ${ONLINE_HEADERS})

target_include_directories(slim_online_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(slim_online_lib PUBLIC slim_common)

# AVX2 support
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(slim_online_lib PRIVATE -mavx2)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(slim_online_lib PRIVATE /arch:AVX2)
endif()

# Executable
add_executable(slim_online src/main.cpp)
target_link_libraries(slim_online PRIVATE slim_online_lib)
