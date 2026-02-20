CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(IJKDUALTABLE)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(SRC_NAME "ijkdualtable")
SET(SRC_DIR "src/${SRC_NAME}")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "-std=c++20")

ADD_EXECUTABLE(testdualtable testdualtable.cpp)

TARGET_INCLUDE_DIRECTORIES(testdualtable PRIVATE "${IJK_DIR}/include")

#---------------------------------------------------------

SET(PROGRAM_FILES ${SRC_DIR}/README
                  ${SRC_DIR}/*.cpp
                  ${SRC_DIR}/*.tpp
                  ${SRC_DIR}/CMakeLists.txt)

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
                      COMMAND tar cvfh ${CMAKE_CURRENT_BINARY_DIR}/ijkdualtable.tar
                                       ${PROGRAM_FILES})

ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
                      COMMAND zip ${CMAKE_CURRENT_BINARY_DIR}/ijkdualtable.zip
                                  ${PROGRAM_FILES})
