project(Google_tests)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# adding the Google_Tests_run target
add_executable(Google_Tests_run
        SimulationSetupTest.cpp
        BoundaryConditionTest.cpp)

target_compile_definitions(Google_Tests_run PRIVATE TEST_BUILD)


target_link_libraries(Google_Tests_run gtest gtest_main)


