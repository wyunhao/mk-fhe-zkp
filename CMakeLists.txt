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
