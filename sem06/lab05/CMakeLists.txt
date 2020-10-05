cmake_minimum_required(VERSION 3.15)
project(lab05 CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
	mainwindow.hpp
	model.hpp
	plot.hpp
)

set(
	SOURCES
	mainwindow.cpp
	model.cpp
	plot.cpp
)

set(
	FORMS
	mainwindow.ui
)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)

add_executable(${PROJECT_NAME} main.cpp ${HEADERS} ${SOURCES} ${FORMS})

target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Widgets)
target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Charts)
