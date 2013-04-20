PROJECT(ViennaFVM)

cmake_minimum_required(VERSION 2.6)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

#specify custom boost path if needed:
#SET(BOOSTPATH "/opt/boost")

# alternate boost path for macs
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   INCLUDE_DIRECTORIES("/opt/local/include")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


#specify include and source directory
INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("../ViennaGrid/")
INCLUDE_DIRECTORIES("../ViennaData/")
INCLUDE_DIRECTORIES("../ViennaMath/")
INCLUDE_DIRECTORIES("../ViennaCL/")
INCLUDE_DIRECTORIES(${BOOSTPATH})

#add definitions, compiler switches, etc.
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
   ADD_DEFINITIONS(-Wall -O0 -g -lpthread)
   #ADD_DEFINITIONS(-Wall -pedantic -O3)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

#list all source files here
ADD_EXECUTABLE(poisson_1d examples/tutorial/poisson_1d.cpp)
ADD_EXECUTABLE(poisson_2d examples/tutorial/poisson_2d.cpp)
ADD_EXECUTABLE(poisson_3d examples/tutorial/poisson_3d.cpp)
ADD_EXECUTABLE(mosfet     examples/tutorial/mosfet.cpp)


##Compatibility with Qt-Creator
file( GLOB_RECURSE QtCreatorCompatibility_SRC
    viennafvm/*.hpp viennafvm/*.h)
add_custom_target( QtCreatorCompatibility SOURCES ${QtCreatorCompatibility_SRC} )
