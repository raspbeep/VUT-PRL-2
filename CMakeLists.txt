#
# @file   CMakeLists.txt
#
# @author David Bayer \n
#         Faculty of Information Technology \n
#         Brno University of Technology \n
#         ibayer@fit.vutbr.cz
#
# @brief  CMakeLists for lab2
#
# @date   14 February 2024, 14:00 (created) \n
#

cmake_minimum_required(VERSION 3.20)

# project name
project(gol LANGUAGES CXX)

# set C++ standard
set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find and setup MPI package
set(MPI_CXX_SKIP_MPICXX ON)

find_package(MPI REQUIRED COMPONENTS CXX)

# compile bcast
add_executable(gol ${CMAKE_CURRENT_SOURCE_DIR}/gol.cpp)
target_link_libraries(gol MPI::MPI_CXX)

# # make pack_assignment target
# add_custom_target(pack_assignment
#   COMMAND ${CMAKE_COMMAND} -E tar "cf" "${CMAKE_CURRENT_BINARY_DIR}/PPP-${PROJECT_NAME}.zip" --format=zip
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/doc.pdf"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/bcast.cpp"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/scatter.cpp")

# # make pack_solution target
# add_custom_target(pack_solution
#   COMMAND ${CMAKE_COMMAND} -E tar "cf" "${CMAKE_CURRENT_BINARY_DIR}/PPP-${PROJECT_NAME}.zip" --format=zip
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/doc.pdf"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/bcast.cpp"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/scatter.cpp"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/bcast_solution.cpp"
#                                        "${CMAKE_CURRENT_SOURCE_DIR}/scatter_solution.cpp")

# if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#   set(CMAKE_CXX_FLAGS                  "-Wall -Wextra")
#   set(CMAKE_CXX_FLAGS_DEBUG            "-O0 -g3")
#   set(CMAKE_CXX_FLAGS_RELEASE          "-Ofast -march=native -mtune=native -fassociative-math -flto=auto")
#   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO   "${CMAKE_CXX_FLAGS_RELEASE} -g1")
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
#   if(WIN32)
#     set(CMAKE_CXX_FLAGS                "/W5")
#     set(CMAKE_CXX_FLAGS_DEBUG          "/Od /debug:full")
#     set(CMAKE_CXX_FLAGS_RELEASE        "/fast /MT /QxHost") # -flto=full is well formed here, see docs!
#     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} /debug:minimal")
#   else()
#     set(CMAKE_CXX_FLAGS                "-Wall -Xlinker -lstdc++")
#     set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -debug full")
#     set(CMAKE_CXX_FLAGS_RELEASE        "-ipo -O3 -no-prec-div -fp-model=fast=2 -xHost")
#     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -debug minimal")
#   endif()
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
#   if(WIN32)
#     set(CMAKE_CXX_FLAGS                "/W5")
#     set(CMAKE_CXX_FLAGS_DEBUG          "/Od /debug:full")
#     set(CMAKE_CXX_FLAGS_RELEASE        "-flto=full /fast /MT /QxHost") # -flto=full is well formed here, see docs!
#     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} /debug:minimal")
#   else()
#     set(CMAKE_CXX_FLAGS                "-Wall -Xlinker -lstdc++")
#     set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -debug full")
#     set(CMAKE_CXX_FLAGS_RELEASE        "-flo=full -fast -xHost")
#     set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -debug minimal")
#   endif()
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
#   set(CMAKE_CXX_FLAGS                  "/W4 /Zc:__cplusplus /EHsc")
#   set(CMAKE_CXX_FLAGS_DEBUG            "/Od /DEBUG:FULL")
#   set(CMAKE_CXX_FLAGS_RELEASE          "/Ox /fp:fast /GL")
#   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO   "${CMAKE_CXX_FLAGS_RELEASE} /DEBUG:FULL")

#   set(CMAKE_SHARED_LINKER_FLAGS        "/LTCG")
# else()
#   message(WARNING "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
# endif()
