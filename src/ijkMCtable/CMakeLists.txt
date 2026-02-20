CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(IJKMCTABLE)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(SRC_NAME "ijkMCtable")
SET(SRC_DIR "src/${SRC_NAME}")
SET(DOXYGEN_CONFIG_FILE "${SRC_NAME}_doxygen.config")
SET(TAR_FILE "${SRC_NAME}.tar")
SET(ZIP_FILE "${SRC_NAME}.zip")
GET_FILENAME_COMPONENT(DEFAULT_INSTALL_PREFIXA "${CMAKE_CURRENT_BINARY_DIR}/.." ABSOLUTE)
GET_FILENAME_COMPONENT(DEFAULT_INSTALL_PREFIXB "${CMAKE_CURRENT_BINARY_DIR}/../.." ABSOLUTE)
SET(DEFAULT_DOC_DIRA "../doc")
SET(DEFAULT_DOC_DIRB "../../doc")
SET(DEFAULT_DOC_DIRC "$ENV{HOME}/doc")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_LIBRARIES(expat)

ADD_COMPILE_OPTIONS(-Wall -Wno-unused-function)

ADD_EXECUTABLE(ijkgenMCtable ijkgenMCtable.cpp ijkgenMCpatch.cpp 
                             ijkMCtable.cpp ijkMCtable_poly.cpp 
                             ijkMCtable_xitIO.cpp clarkson_hull.c)
ADD_EXECUTABLE(ijkMCtable_invert ijkMCtable_invert.cpp
                                 ijkMCtable.cpp ijkMCtable_poly.cpp
                                 ijkMCtable_xitIO.cpp)
ADD_EXECUTABLE(ijkdiffMCtable ijkdiffMCtable.cpp ijkMCtable.cpp 
                              ijkMCtable_poly.cpp ijkMCtable_xitIO.cpp)
ADD_EXECUTABLE(ijkMCtableinfo ijkMCtableinfo.cpp ijkMCtable.cpp 
                              ijkMCtable_poly.cpp ijkMCtable_xitIO.cpp)
ADD_EXECUTABLE(outambigtable outambigtable.cpp 
                             ijkMCtable_poly.cpp ijkMCtable_ambig.cpp)

# CURRENTLY NOT USED
# ADD_EXECUTABLE(ijktext2c ijktext2c.cpp)


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
   INSTALL(TARGETS ijkgenMCtable DESTINATION "bin")
   INSTALL(TARGETS ijkMCtable_invert DESTINATION "bin")
   INSTALL(TARGETS ijkdiffMCtable DESTINATION "bin")
   INSTALL(TARGETS ijkMCtableinfo DESTINATION "bin")
   INSTALL(TARGETS outambigtable DESTINATION "bin")
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
                  ${SRC_DIR}/*.tpp
                  ${SRC_DIR}/*.c
                  ${SRC_DIR}/*.xsd
                  ${SRC_DIR}/CMakeLists.txt
                  ${SRC_DIR}/${DOXYGEN_CONFIG_FILE})

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
  COMMAND tar cvfh ${CMAKE_CURRENT_BINARY_DIR}/${TAR_FILE}
                   ${PROGRAM_FILES})

ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
  COMMAND zip ${CMAKE_CURRENT_BINARY_DIR}/${ZIP_FILE}
              ${PROGRAM_FILES})

#---------------------------------------------------------

if (NOT DEFINED DOC_DIR)
  if (IS_DIRECTORY "${DEFAULT_DOC_DIRA}")
    SET(DOC_DIR "${DEFAULT_DOC_DIRA}")
  ELSEIF (IS_DIRECTORY "${DEFAULT_DOC_DIRB}")
    SET(DOC_DIR "${DEFAULT_DOC_DIRB}")
  ELSEIF (IS_DIRECTORY "${DEFAUL_DOC_DIRC}")
    SET(DOC_DIR "${DEFAULT_DOC_DIRC}")
  ELSE()
    MESSAGE(WARNING " Warning: Unable to determine directory for doxygen files.\n"
                    " Use -DDOC_DIR=... to specify directory containing doxygen documentation.\n"
                    " Create directory ${DEFAULT_DOC_DIRA} or ${DEFAULT_DOC_DIRB}"
                    " and rerun cmake to enable doxygen if you want cmake"
                    " to automatically use one of those directories as the"
                    " doxygen directory.\n"
                    " Target doc disabled.")
  ENDIF()
ENDIF()

IF (DEFINED DOC_DIR)
  IF (IS_DIRECTORY "${DOC_DIR}")
    GET_FILENAME_COMPONENT(DOC_DIR_ABSOLUTE_PATH
        "${CMAKE_CURRENT_BINARY_DIR}/${DOC_DIR}" ABSOLUTE)
    SET(DOXYGEN_GENERATED_FILE "${SRC_NAME}_doxygen.generated")
    FILE(READ "${DOXYGEN_CONFIG_FILE}" DOXYGEN_CONFIG)
    FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${DOXYGEN_GENERATED_FILE}
               "${DOXYGEN_CONFIG}\nINPUT=${IJK_DIR}/${SRC_DIR}\nOUTPUT_DIRECTORY=${DOC_DIR_ABSOLUTE_PATH}")
    ADD_CUSTOM_TARGET(doc WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                          COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/${DOXYGEN_GENERATED_FILE})
    MESSAGE(STATUS "Directory containing doxygen documentation: ${DOC_DIR_ABSOLUTE_PATH}")
  ELSE()
    # DOC_DIR was defined in command line but does not exist.
    MESSAGE(WARNING "\"${DOC_DIR}\" does not exist or is not a directory.\n"
                    " Directory \"DIR_NAME\" must exist to use -DDOC_DIR=\{DIR\_NAME\}.\n"
                    " Target doc disabled.")
  ENDIF()
ENDIF()

