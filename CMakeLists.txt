cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(AlignmentRepaC)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(../AlignmentC ../AlignmentC_build)

include_directories(../rapidjson/include ../AlignmentC)

add_library(libAlignmentRepaC AlignmentRepa.cpp AlignmentAesonRepa.cpp AlignmentRandomRepa.cpp AlignmentPracticableRepa.cpp AlignmentPracticableIORepa.cpp)

add_executable(AlignmentRepaC main.cpp)

target_link_libraries(AlignmentRepaC PUBLIC libAlignmentC libAlignmentRepaC)

target_include_directories(AlignmentRepaC PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          )

install(TARGETS AlignmentRepaC DESTINATION bin)
install(TARGETS libAlignmentRepaC DESTINATION lib)
install(FILES AlignmentRepa.h AlignmentAesonRepa.h AlignmentRandomRepa.h AlignmentPracticableRepa.h AlignmentPracticableIORepa.h
  DESTINATION include
  )
