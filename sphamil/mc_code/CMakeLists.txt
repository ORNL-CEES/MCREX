CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Get Trilinos as one entity
FIND_PACKAGE(
        Trilinos 
        PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH} 
        REQUIRED)

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("End of Trilinos details\n")

MESSAGE("Trilinos Flags")
MESSAGE(${Trilinos_CXX_COMPILER_FLAGS})

SET(CMAKE_C_COMPILER        ${Trilinos_C_COMPILER})
SET(CMAKE_CXX_COMPILER      ${Trilinos_CXX_COMPILER})
SET(CMAKE_Fortran_COMPILER  ${Trilinos_Fortran_COMPILER})

PROJECT(MCREX)

SET(SOURCES
  AdditiveSchwarzWrapper.cc
  BelosSolver.cc
  ChebyshevIteration.cc
  ChebyshevPolynomial.cc
  GmresPolynomial.cc
  LinearSolverFactory.cc
  LinearSystem.cc
  LinearSystemFactory.cc
  MCREX_Solver.cc
  MC_Data.cc
  MonteCarloSolver.cc
  NeumannPolynomial.cc
  Polynomial.cc
  PolynomialBasis.cc
  PolynomialFactory.cc
  PolynomialPreconditioner.cc
  PolynomialUtils.cc
  RichardsonIteration.cc
  SyntheticAcceleration.cc
  )

FILE(GLOB HEADERS *.hh)

# Enable most warnings, explicitly disable:
#  unused-parameter: Belos results in many of these
#  long-long: gtest results in many of these
#  no-variadic-macros: gtest again
SET(LOCAL_CXX_FLAGS 
    "-ansi -pedantic -Wall -Wextra -Wno-unused-parameter -Wno-long-long -Wno-variadic-macros"
   )

SET(CXX11_FLAGS "")
IF( CMAKE_COMPILER_IS_GNUCC )
    MESSAGE("CMAKE_CXX_COMPILER_VERSION")
    MESSAGE(${CMAKE_CXX_COMPILER_VERSION})
    IF( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7.0 )
        SET(CXX11_FLAGS "-std=c++0x")
    ELSE()
        SET(CXX11_FLAGS "-std=c++11")
    ENDIF()
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS} ${LOCAL_CXX_FLAGS} ${CXX11_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trilinos_C_COMPILER_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS}")
MESSAGE("CMAKE_CXX_FLAGS")
MESSAGE(${CMAKE_CXX_FLAGS})

INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

ADD_LIBRARY(mcrex ${SOURCES})

ADD_EXECUTABLE(main mcrex.cc)
ADD_EXECUTABLE(multirun multirun.cc)

TARGET_LINK_LIBRARIES(main     mcrex ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 
TARGET_LINK_LIBRARIES(multirun mcrex ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 

# Testing
ENABLE_TESTING()

ADD_SUBDIRECTORY(gtest)
INCLUDE_DIRECTORIES(gtest/include)

#CONFIGURE_FILE(testMonteCarloSolver.xml ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

SET(TESTS
  testMonteCarloSolver
  testChebyshevIteration
  testLinearSolverFactory
  testLinearSystemFactory
  testMC_Data
  testMCSA
  testPolynomial
  testPolynomialBasis
  testPolynomialUtils
  testRichardsonIteration
  )

FOREACH( testname ${TESTS} )
    ADD_EXECUTABLE(${testname} 
            ${testname}.cc )
    TARGET_LINK_LIBRARIES(${testname} 
            mcrex 
            gtest_all
            gtest_main 
            ${Trilinos_LIBRARIES} 
            ${Trilinos_TPL_LIBRARIES}) 
    ADD_TEST(${testname} 
            ${CMAKE_CURRENT_BINARY_DIR}/${testname} )
ENDFOREACH()

SET(PARALLEL_TESTS
    testLinearSystemFactory
   )

FOREACH( testname ${PARALLEL_TESTS} )
    ADD_TEST(NAME ${testname}_2 COMMAND mpirun -np 2 ${CMAKE_CURRENT_BINARY_DIR}/${testname} )
    ADD_TEST(NAME ${testname}_4 COMMAND mpirun -np 4 ${CMAKE_CURRENT_BINARY_DIR}/${testname} )
ENDFOREACH()

