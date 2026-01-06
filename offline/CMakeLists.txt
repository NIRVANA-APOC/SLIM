# Offline index builder sources
set(OFFLINE_SOURCES
    src/builder.cpp
    src/hierarchical.cpp
    src/encoder.cpp
)

set(OFFLINE_HEADERS
    include/slim/offline/builder.h
    include/slim/offline/hierarchical.h
    include/slim/offline/encoder.h
)

add_library(slim_offline_lib STATIC ${OFFLINE_SOURCES} ${OFFLINE_HEADERS})

target_include_directories(slim_offline_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(slim_offline_lib PUBLIC slim_common)

# Executable
add_executable(slim_offline src/main.cpp)
target_link_libraries(slim_offline PRIVATE slim_offline_lib)
