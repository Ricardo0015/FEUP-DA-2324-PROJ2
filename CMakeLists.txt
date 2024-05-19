cmake_minimum_required(VERSION 3.28)
project(FEUP_DA_2324_PROJ2)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(FEUP_DA_2324_PROJ2
        Graph.h
        Heap.cpp
        Heap.h
        main.cpp
        MutablePriorityQueue.h
        UFDS.cpp
        UFDS.h
        Management.h
        Management.cpp
        points.cpp
        points.h
)
