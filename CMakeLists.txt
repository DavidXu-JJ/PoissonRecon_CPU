cmake_minimum_required(VERSION 3.17)
project(Practice)

set(CMAKE_CXX_STANDARD 14)

add_executable(Practice main.cpp Polynomial.h Factor.cpp Factor.h PPolynomial.h FunctionData.h BinaryNode.h SparseMatrix.h Vector.h Allocator.h Geometry.h Geometry.cpp MarchingCubes.cpp MarchingCubes.h Octree.h)