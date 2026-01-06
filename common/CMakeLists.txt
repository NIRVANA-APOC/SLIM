set(COMMON_SOURCES
    src/graph.cpp
    src/index.cpp
    src/spectral/laplacian.cpp
    src/spectral/eigen.cpp
)

set(COMMON_HEADERS
    include/slim/graph.h
    include/slim/index.h
    include/slim/spectral/laplacian.h
    include/slim/spectral/eigen.h
)

add_library(slim_common STATIC ${COMMON_SOURCES} ${COMMON_HEADERS})

target_include_directories(slim_common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(slim_common PUBLIC cxx_std_17)

if(OpenMP_CXX_FOUND)
    target_link_libraries(slim_common PUBLIC OpenMP::OpenMP_CXX)
endif()

# Optional allocators
if(TCMALLOC_LIB)
    target_link_libraries(slim_common PUBLIC ${TCMALLOC_LIB})
elseif(JEMALLOC_LIB)
    target_link_libraries(slim_common PUBLIC ${JEMALLOC_LIB})
endif()
