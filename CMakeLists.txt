cmake_minimum_required(VERSION 3.10)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${project_SOURCE_DIR}/cmake")
project(prp)

include(ExternalProject)
find_package(Git REQUIRED)


ExternalProject_Add(
    args
    PREFIX ${CMAKE_BINARY_DIR}/args
    GIT_REPOSITORY https://github.com/Taywee/args
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
)
ExternalProject_Get_Property(args source_dir)
set(ARGS_INCLUDE_DIR ${source_dir} CACHE INTERNAL "Path to include folder for args")

find_package(SDL2 REQUIRED)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall --pedantic -O3 -g -finline-functions -march=native")
#set(CMAKE_CXX_FLAGS "-Wall --pedantic -g -finline-functions -march=native")

set(LIB_SOURCE_FILES)

add_executable(prp src/main.cc ${LIB_SOURCE_FILES})
add_dependencies(prp args)
target_include_directories(prp PUBLIC ${ARGS_INCLUDE_DIR} ${SDL2_INCLUDE_DIRS})
target_link_libraries(prp ${SDL2_LIBRARIES})


# TESTING
ExternalProject_Add(
    doctest
    PREFIX ${CMAKE_BINARY_DIR}/doctest
    GIT_REPOSITORY https://github.com/onqtam/doctest.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
)
ExternalProject_Get_Property(doctest source_dir)
set(DOCTEST_INCLUDE_DIR ${source_dir}/doctest CACHE INTERNAL "Path to include folder for doctest")
enable_testing()

add_executable(prp_test tests/prp.test.cc ${LIB_SOURCE_FILES} tests/doctest.cc)
add_dependencies(prp_test doctest)
target_include_directories(prp_test PUBLIC ${DOCTEST_INCLUDE_DIR})
add_test(prp prp_test)
