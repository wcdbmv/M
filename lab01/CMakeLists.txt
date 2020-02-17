cmake_minimum_required(VERSION 3.15)
project(lab01)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(
	-Werror

	-Wall
	-Wextra
	-Wpedantic

	-Wcast-align
	-Wcast-qual
	-Wconversion
	-Wctor-dtor-privacy
	-Wenum-compare
	-Wfloat-equal
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Woverloaded-virtual
	-Wredundant-decls
	-Wsign-conversion
	-Wsign-promo
)

set(
	HEADERS
	polynomial/polynomial.hpp)

set(
	SOURCES
	polynomial/polynomial.cpp)

add_executable(lab01 main.cpp ${HEADERS} ${SOURCES})
