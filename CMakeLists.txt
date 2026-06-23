cmake_minimum_required(VERSION 3.8)
project(PRISM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Directory for parameters
set(PARAMS_DIR "${CMAKE_SOURCE_DIR}/params")  # or CMAKE_INSTALL_PREFIX based path

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto")

SET(SOURCE
  src/W_final.cc 
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
# configure_file(RNA.pc.in RNA.pc @ONLY)

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
