cmake_minimum_required(VERSION 3.15)
project(lab03 CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5 COMPONENTS PrintSupport REQUIRED)

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

include_directories(../qcustomplot)

set(
	HEADERS
	../qcustomplot/qcustomplot.h
	mainwindow.hpp
	solve.hpp
)

set(
	SOURCES
	../qcustomplot/qcustomplot.cpp
	mainwindow.cpp
	solve.cpp
)

set(
	FORMS
	mainwindow.ui
)

add_executable(${PROJECT_NAME} main.cpp ${HEADERS} ${SOURCES} ${FORMS})

target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Widgets)
target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::PrintSupport)
