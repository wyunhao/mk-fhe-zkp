cmake_minimum_required(VERSION 3.13)
project(mk_rlwe_prototype LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(SEAL 4.1 REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/seal-install/lib/cmake/SEAL-4.1" NO_DEFAULT_PATH)

add_executable(prototype prototype.cpp)
target_link_libraries(prototype PRIVATE SEAL::seal)
target_compile_options(prototype PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(seal_bench seal_bench.cpp)
target_link_libraries(seal_bench PRIVATE SEAL::seal)
target_compile_options(seal_bench PRIVATE -O2 -Wall)

add_executable(seal_bench_bgv_bfv seal_bench_bgv_bfv.cpp)
target_link_libraries(seal_bench_bgv_bfv PRIVATE SEAL::seal)
target_compile_options(seal_bench_bgv_bfv PRIVATE -O2 -Wall)

add_executable(microbench microbench.cpp)
target_link_libraries(microbench PRIVATE SEAL::seal)
target_compile_options(microbench PRIVATE -O2 -Wall)

add_executable(prototype_ext prototype_ext.cpp)
target_link_libraries(prototype_ext PRIVATE SEAL::seal)
target_compile_options(prototype_ext PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(prototype_ptct prototype_ptct.cpp)
target_link_libraries(prototype_ptct PRIVATE SEAL::seal)
target_compile_options(prototype_ptct PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(prototype_full_chain prototype_full_chain.cpp)
target_link_libraries(prototype_full_chain PRIVATE SEAL::seal)
target_compile_options(prototype_full_chain PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(prototype_full_chain_nozk prototype_full_chain_nozk.cpp)
target_link_libraries(prototype_full_chain_nozk PRIVATE SEAL::seal)
target_compile_options(prototype_full_chain_nozk PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(prototype_longest_chain prototype_longest_chain.cpp)
target_link_libraries(prototype_longest_chain PRIVATE SEAL::seal)
target_compile_options(prototype_longest_chain PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(prototype_tight_chain prototype_tight_chain.cpp)
target_link_libraries(prototype_tight_chain PRIVATE SEAL::seal)
target_compile_options(prototype_tight_chain PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(prototype_longest_chain_p17 prototype_longest_chain_p17.cpp)
target_link_libraries(prototype_longest_chain_p17 PRIVATE SEAL::seal)
target_compile_options(prototype_longest_chain_p17 PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

add_executable(proto_large_pl_tighten proto_large_pl_tighten.cpp)
target_link_libraries(proto_large_pl_tighten PRIVATE SEAL::seal)
target_compile_options(proto_large_pl_tighten PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)

# util library (reusable helpers)
add_library(mk_util STATIC util/ring_switch.cpp)
target_include_directories(mk_util PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mk_util PUBLIC SEAL::seal)
target_compile_options(mk_util PRIVATE -O2 -Wall)

add_executable(prototype_ring_switch prototype_ring_switch.cpp)
target_link_libraries(prototype_ring_switch PRIVATE mk_util SEAL::seal)
target_compile_options(prototype_ring_switch PRIVATE -O2 -Wall -Wno-unused-function -Wno-unused-variable)
