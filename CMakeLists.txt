
# mkdir build
# cd build
# cmake ..
# cmake --build . --parallel
# ctest --output-on-failure

# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10) 

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(SPINAS)
include(CTest)
include(GNUInstallDirs)

# Set precision type
set(PRECISION_TYPE "LONG_DOUBLE" CACHE STRING "The precision type to be used: LONG_DOUBLE")


# Set compiler flags
if(${PRECISION_TYPE} STREQUAL "LONG_DOUBLE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DWITH_LONG_DOUBLE -O3 -fPIC")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC")
endif()

# Set directories
set(SRC_DIR source)
set(INC_DIR include)
set(OBJ_DIR obj)
set(LIB_DIR lib) 

# Include directories
include_directories(${INC_DIR})

# Add main library
file(GLOB_RECURSE LIB_SRC_FILES "${SRC_DIR}/*.cpp")
add_library(spinas_shared SHARED ${LIB_SRC_FILES})
set_target_properties(spinas_shared PROPERTIES OUTPUT_NAME spinas_shared)
add_library(spinas_static STATIC ${LIB_SRC_FILES}) 
set_target_properties(spinas_static PROPERTIES OUTPUT_NAME spinas)


####################################################
#                      Core Tests
####################################################
if(BUILD_TESTING)
    file(GLOB TEST_SRCS tests/*.cpp)
    add_executable(test_spinas ${TEST_SRCS})
    target_link_libraries(test_spinas spinas_static)
    add_test(NAME test_spinas COMMAND test_spinas)
endif()

####################################################
#                      SM Tests
####################################################
set(SM_DIR SM)
file(GLOB_RECURSE SM_SRC_FILES "${SM_DIR}/*.cpp")
list(REMOVE_ITEM SM_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/SM/test.cpp")
add_library(spinas_SM_Obj OBJECT ${SM_SRC_FILES})  # create object library
add_library(spinas_SM SHARED $<TARGET_OBJECTS:spinas_SM_Obj>)  # create shared library from object files
target_link_libraries(spinas_SM spinas_static)  # link the shared library
if(BUILD_TESTING)
    add_executable(test_SM SM/test.cpp tests/test_utils.cpp $<TARGET_OBJECTS:spinas_SM_Obj>)
    target_include_directories(test_SM PRIVATE tests)
    target_link_libraries(test_SM spinas_static)
    add_test(NAME test_SM COMMAND test_SM)
endif()


####################################################
#                      Installation
####################################################
# If installing in a non-default location, use e.g.: cmake -DCMAKE_INSTALL_PREFIX=${HOME}/.local ..
install(TARGETS spinas_shared spinas_static
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
#install(DIRECTORY ${INC_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
