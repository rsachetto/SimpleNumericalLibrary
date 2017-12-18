list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

FIND_PACKAGE(Criterion)


if(CRITERION_FOUND)
    set(SOURCE_FILES_TEST tests.cpp)
    add_executable(TestSuite ${SOURCE_FILES_TEST})
    ADD_DEPENDENCIES(TestSuite simple_numerical_library)
    target_link_libraries(TestSuite simple_numerical_library m criterion)
endif()