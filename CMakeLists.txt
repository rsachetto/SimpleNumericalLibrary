cmake_minimum_required(VERSION 3.0.1)
project(SimpleNumericalLibrary)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)



SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )


add_library(simple_numerical_library SHARED
        src/DirectSolver.cpp
        src/DirectSolver.h
        src/Function.h
        src/GaussQuadrature.cpp
        src/GaussQuadrature.h
        src/IdentityMatrix.cpp
        src/IdentityMatrix.h
        src/IterativeSolver.cpp
        src/IterativeSolver.h
        src/Matrix.cpp
        src/Matrix.h
        src/MatrixFileHandler.cpp
        src/MatrixFileHandler.h
        src/NonLinearSolver.cpp
        src/NonLinearSolver.h
        src/RootFinder.cpp
        src/RootFinder.h
        src/TridiagonalMatrix.cpp
        src/TridiagonalMatrix.h
        src/TridiagonalSolver.cpp
        src/TridiagonalSolver.h
        src/util.h
        src/Vector.cpp
        src/Vector.h
        src/VectorFileHandler.cpp
        src/VectorFileHandler.h
        )

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

ADD_SUBDIRECTORY(tests)