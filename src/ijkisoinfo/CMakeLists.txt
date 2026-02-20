CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(IJKISOINFO)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(SRC_NAME "ijkisoinfo")
SET(SRC_DIR "src/${SRC_NAME}")
SET(LIB_DIR_LIST "${CMAKE_CURRENT_BINARY_DIR}/../lib"
                 "${CMAKE_CURRENT_BINARY_DIR}/../../lib"
                 "${IJK_DIR}/lib"
                 "$ENV{HOME}/lib"
                 "/opt/lib")
GET_FILENAME_COMPONENT(DEFAULT_INSTALL_PREFIXA "${CMAKE_CURRENT_BINARY_DIR}/.." ABSOLUTE)
GET_FILENAME_COMPONENT(DEFAULT_INSTALL_PREFIXB "${CMAKE_CURRENT_BINARY_DIR}/../.." ABSOLUTE)

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")

# Possibly turn on usage of zlib compression (requires linking with libz)
# (i.e., programs compiled with ITKNrrdIO must also be compiled with zlib)
ADD_DEFINITIONS(-DTEEM_ZLIB=1)

# Turn on TEEM_BUILD so that the proper dll export def's are 
# used on windows builds.
ADD_DEFINITIONS(-DTEEM_BUILD=1)

SET(CMAKE_CXX_FLAGS "-std=c++20")

ADD_EXECUTABLE(ijkisoinfo ijkisoinfo.cpp)

#---------------------------------------------------------

FIND_LIBRARY(NRRDIO_LIB NAMES NrrdIO HINTS ${LIB_DIR_LIST})

IF (NOT NRRDIO_LIB)
   MESSAGE(FATAL_ERROR "Unable to find library NrrdIO.\n"
                       "Source code and CMakeLists.txt to create library NrrdIO is in ${IJK_DIR}/src/NrrdIO.\n"
                       "Library NrrdIO should be in one of the following directories:\n"
                       "${LIB_DIR_LIST}")
ELSE()
   MESSAGE(STATUS "Found NrrdIO library: ${NRRDIO_LIB}.")
   TARGET_LINK_DIRECTORIES(ijkisoinfo PRIVATE ${LIB_DIR_LIST})
   TARGET_LINK_LIBRARIES(ijkisoinfo NrrdIO z)
ENDIF()

#---------------------------------------------------------

IF (NOT DEFINED IJK_INSTALL_PREFIX)
   IF (EXISTS "${DEFAULT_INSTALL_PREFIXA}/bin")
      SET(IJK_INSTALL_PREFIX "${DEFAULT_INSTALL_PREFIXA}")
   ELSEIF (EXISTS "${DEFAULT_INSTALL_PREFIXB}/bin")
      SET(IJK_INSTALL_PREFIX "${DEFAULT_INSTALL_PREFIXB}")
   ELSEIF (EXISTS "${IJK_DIR}/bin")
      SET(IJK_INSTALL_PREFIX "${IJK_DIR}")
   ENDIF()
ENDIF()

IF (DEFINED IJK_INSTALL_PREFIX)
   SET(CMAKE_INSTALL_PREFIX "${IJK_INSTALL_PREFIX}")
   MESSAGE(STATUS "Install directory: ${IJK_INSTALL_PREFIX}/bin")
   INSTALL(TARGETS ijkisoinfo DESTINATION "bin")
ELSE()
   MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                   " Use -DIJK_INSTALL_PREFIX=... to specify install directory.\n"
                   " Create directory ${DEFAULT_INSTALL_PREFIXA}/bin or"
                   " ${DEFAULT_INSTALL_PREFIXB}/bin and rerun cmake"
                   " if you want cmake to automatically use one of those directories"
                   " as the install directory.")
ENDIF()

#---------------------------------------------------------

SET(PROGRAM_FILES ${SRC_DIR}/README
                  ${SRC_DIR}/*.cpp
                  ${SRC_DIR}/*.h
                  ${SRC_DIR}/CMakeLists.txt)

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
                      COMMAND tar cvfh ${CMAKE_CURRENT_BINARY_DIR}/ijkisoinfo.tar
                      ${PROGRAM_FILES})

ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
                      COMMAND zip ${CMAKE_CURRENT_BINARY_DIR}/ijkisoinfo.zip
                      ${PROGRAM_FILES})
