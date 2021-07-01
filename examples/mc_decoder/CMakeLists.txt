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
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/codec_polar.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Decoder_polar_SCL_mcfast_sys.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Decoder_polar_SCL_mcfast_sys.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Decoder_polar_SCL_oldfast_sys.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Decoder_polar_SCL_oldfast_sys.cpp
                    )
# Add glog to project
find_package (glog 0.5.0 REQUIRED)
target_link_libraries (main PUBLIC glog::glog)

# Link with the "Threads library (required to link with AFF3CT after)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Link with AFF3CT
set (AFF3CT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
find_package(AFF3CT CONFIG 2.3.5 REQUIRED)
target_link_libraries(main PRIVATE aff3ct::aff3ct-static-lib)