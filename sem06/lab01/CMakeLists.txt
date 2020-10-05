cmake_minimum_required(VERSION 3.15)
project(lab01)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

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
	cauchy/cauchy.hpp
	polynomial/polynomial.hpp
	range/range.hpp
)

set(
	SOURCES
	cauchy/cauchy.cpp
	polynomial/polynomial.cpp
	range/range.cpp
)

add_executable(lab01 main.cpp ${HEADERS} ${SOURCES})
