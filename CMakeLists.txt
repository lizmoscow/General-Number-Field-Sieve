cmake_minimum_required(VERSION 3.6)
project(General_Number_Field_Sieve)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES "main.cpp"
        "Initialization.cpp"
        "Sieving.cpp"
        "LinearDependency.cpp"
        "Chineese remainder.cpp"
        "GNFS.cpp")

add_executable(General_Number_Field_Sieve ${SOURCE_FILES})

target_link_libraries(General_Number_Field_Sieve ntl gmp armadillo givaro linbox)