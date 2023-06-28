# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)
FetchContent_Declare(
    _scale_straw
    GIT_REPOSITORY https://github.com/aidenlab/straw.git
    GIT_TAG        2525edc29bbb48463799cad94cbd6e5e810210a0
)

FetchContent_MakeAvailable(_scale_straw)

find_package(ZLIB REQUIRED)
find_package(CURL REQUIRED)

add_library(_scale_strawlib)
target_sources(_scale_strawlib PRIVATE ${_scale_straw_SOURCE_DIR}/C++/straw.cpp)
target_include_directories(_scale_strawlib PUBLIC ${_scale_straw_SOURCE_DIR}/C++/)
target_link_libraries(_scale_strawlib PUBLIC CURL::libcurl ZLIB::ZLIB)

add_library(scale)

target_sources(scale PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/d_ppNormVector.cpp ${CMAKE_CURRENT_SOURCE_DIR}/d_thMul.cpp ${CMAKE_CURRENT_SOURCE_DIR}/finito.cpp)

target_include_directories(scale PUBLIC include/)


add_executable(scale_tools)
set_target_properties(scale_tools PROPERTIES OUTPUT_NAME scale)

target_sources(scale_tools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/finitoMain.cpp ${CMAKE_CURRENT_SOURCE_DIR}/getSingleMatrix.cpp ${CMAKE_CURRENT_SOURCE_DIR}/getGWMatrix.cpp)
target_link_libraries(scale_tools PUBLIC scale _scale_strawlib)


include(GNUInstallDirs)
install(
        TARGETS scale scale_tools
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

