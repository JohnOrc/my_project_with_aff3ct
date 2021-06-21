cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0054 NEW)

project (main)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify bin path
set (EXECUTABLE_OUTPUT_PATH bin/)

# Create the executable from sources
add_executable(main ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Codec_polar_mc.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Codec_polar_mc.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Decoder_polar_mc.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Decoder_polar_mc.cpp
                    )

# Link with the "Threads library (required to link with AFF3CT after)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Link with AFF3CT
set (AFF3CT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
find_package(AFF3CT CONFIG 2.3.5 REQUIRED)
target_link_libraries(main PRIVATE aff3ct::aff3ct-static-lib)
