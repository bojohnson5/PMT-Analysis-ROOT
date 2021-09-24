cmake_minimum_required(VERSION 2.8)
project(analysis)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Release)
set(ROOT_DIR "/N/u/coherent/Carbonate/share/root-6.25_build")

find_package(ROOT 6.25 CONFIG REQUIRED)

add_executable(analysis main.cpp rooter.cpp rooter.hpp)

target_link_libraries(analysis PUBLIC ROOT::Core ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Graf3d ROOT::RIO)
