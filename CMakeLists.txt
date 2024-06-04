cmake_minimum_required(VERSION 3.28)
project(MatheProjekt)

set(CMAKE_CXX_STANDARD 26)

add_executable(MatheProjekt main.cpp
        CMyVektor.cpp
        CMyVektor.h
        CMyMatrix.cpp
        CMyMatrix.h
        C_DGLSolver.cpp
        C_DGLSolver.h)
