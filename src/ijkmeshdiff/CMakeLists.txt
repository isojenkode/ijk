CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(IJKMESHDIFF)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(IJKMESHDIFF_DIR "src/ijkmeshdiff")

SET(SRC_NAME "ijkmeshdiff")
SET(SRC_DIR "src/${SRC_NAME}")
SET(TAR_FILE "${SRC_NAME}.tar")
SET(ZIP_FILE "${SRC_NAME}.zip")
GET_FILENAME_COMPONENT(DEFAULT_INSTALL_PREFIXA "${CMAKE_CURRENT_BINARY_DIR}/.." ABSOLUTE)
GET_FILENAME_COMPONENT(DEFAULT_INSTALL_PREFIXB "${CMAKE_CURRENT_BINARY_DIR}/../.." ABSOLUTE)

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")

ADD_COMPILE_OPTIONS(-Wall -Wno-unused-function)

SET(CMAKE_CXX_FLAGS "-std=c++20")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

ADD_EXECUTABLE(ijkmeshdiff ijkmeshdiff.cpp)

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
   INSTALL(TARGETS ijkmeshdiff DESTINATION "bin")
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
                  ${SRC_DIR}/CMakeLists.txt)

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
                      COMMAND tar cvfh ${CMAKE_CURRENT_BINARY_DIR}/${TAR_FILE}
                                       ${PROGRAM_FILES})

ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
                      COMMAND zip ${CMAKE_CURRENT_BINARY_DIR}/${ZIP_FILE}
                                  ${PROGRAM_FILES})
