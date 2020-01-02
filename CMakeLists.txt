cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(AlignmentRepaC)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Threads)

add_subdirectory(../AlignmentC ../AlignmentC_build)

include_directories(../rapidjson/include ../AlignmentC)

add_library(AlignmentRepaC AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp)

add_executable(AlignmentRepaC_test main.cpp)

target_link_libraries(AlignmentRepaC_test PUBLIC AlignmentRepaC AlignmentC Threads::Threads)

target_include_directories(AlignmentRepaC_test PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          )

install(TARGETS AlignmentRepaC DESTINATION lib)
install(FILES AlignmentRepa.h AlignmentAesonRepa.h AlignmentRandomRepa.h AlignmentPracticableRepa.h AlignmentPracticableIORepa.h
  DESTINATION include
  )
