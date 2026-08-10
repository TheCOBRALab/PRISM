cmake_minimum_required(VERSION 3.15)
project(PRISM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Directory for parameters
set(PARAMS_DIR "${CMAKE_SOURCE_DIR}/params")  # or CMAKE_INSTALL_PREFIX based path

# Default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------------------------------------------
#                             C++ Flags
# -------------------------------------------------------------------

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Strict")
    message(STATUS "Building in Debug mode")
    add_compile_options(
        -g
        -O0
        -Wall
        -Wextra
        # -Wpedantic
        # -Wshadow
        # -Wsign-compare
        # -Wconversion
        # -Wsign-conversion
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
        -fno-omit-frame-pointer
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Strict")
        add_compile_options(-Werror)
    endif()

    add_link_options(
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
else()
    # Normal (Release or other non-Debug builds)
    message(STATUS "Building in Release/Normal mode")
    add_compile_options(
        -O3
        -march=native
        -Wall
        -Wextra
    )
endif()

set(SOURCE
  src/pseudo_loop_trace.cc
  src/pseudo_loop.cc 
  src/part_func.cc
  src/PRISM.cc
  src/cmdline.cc
  src/Result.cc	
  src/s_energy_matrix.cc	
  src/Hotspot.cc
  src/sparse_tree.cc
  src/dot_plot.cc
  src/mea.cc
  src/centroid.cc
  src/SHAPE.cc
)

set(vienna_SOURCE
    src/ViennaRNA/params/default.cc
    src/ViennaRNA/params/params.cc
    src/ViennaRNA/params/io.cc
    src/ViennaRNA/model.cc 
    src/ViennaRNA/utils.cc
)

add_library(RNA STATIC ${vienna_SOURCE})

target_include_directories(RNA PRIVATE .)

add_executable(PRISM ${SOURCE})

# Pass it as a compile definition
target_compile_definitions(PRISM PRIVATE 
    PARAMS_DIR="${PARAMS_DIR}"
)

# Install CParty into the bin directory
install(TARGETS PRISM RUNTIME DESTINATION bin)

target_link_libraries(PRISM PRIVATE RNA)

include_directories(src)
