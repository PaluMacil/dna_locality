cmake_minimum_required(VERSION 3.17)
project(dna_locality)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
add_executable(${PROJECT_NAME} src/main.cpp src/SimilarityMatrix.cpp src/SimilarityMatrix.h src/LocatedScore.cpp src/LocatedScore.h)