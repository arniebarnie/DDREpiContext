cmake_minimum_required(VERSION 3.10)
project(DDREpiContext)

set(CMAKE_CXX_STANDARD 17)

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})

find_package(Boost COMPONENTS REQUIRED filesystem iostreams system regex python)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_executable(DDREpiContext main.cpp util.h hmm.h)

TARGET_LINK_LIBRARIES(DDREpiContext LINK_PUBLIC ${Boost_LIBRARIES})