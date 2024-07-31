cmake_minimum_required(VERSION 3.28)
project(MDS)

#set(CMAKE_CXX_STANDARD 23)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# Include ZeroMQ
# set(ZMQ_INCLUDE_DIR "C:/msys64/mingw64/include")
# set(ZMQ_LIB_DIR "C:/msys64/mingw64/lib")
# include_directories(${ZMQ_INCLUDE_DIR})
# link_directories(${ZMQ_LIB_DIR})

set(EIGEN_PATH "C:/eigen-3.4.0")
include_directories(${EIGEN_PATH})

include_directories(${PROJECT_SOURCE_DIR}/Headers)
include_directories(.)

# Define the ABI macro for the whole project
# add_compile_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

set(MATPLOTPP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(MATPLOTPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory("External Non-compiled Libraries/matplotplusplus")



# inih library
add_library(inih STATIC "External Non-compiled Libraries/inih/ini.c")
# target_include_directories(inih PUBLIC third_party/inih)
# target_compile_features(inih PUBLIC cxx_std_23)

# Source files
set(SOURCES
        abstract_factory.cpp
        boundary_conditions.cpp
        builder.cpp
        command.cpp
        concretecommands.cpp
        datasaver.cpp
        director.cpp
        errorcodes.cpp
        leapfrog.cpp
        main.cpp
        menu.cpp
        progressbar.cpp
        simulation.cpp
        simulation2dap.cpp
        simulation3dcs.cpp
        Headers/core/constants.h
        Headers/core/errorcodes.h
        Headers/core/types.h
        Headers/core/concretecommands.h
        Headers/core/command.h
        Headers/simulation/abstract_factory.h
        Headers/simulation/boundary_conditions.h
        Headers/simulation/builder.h
        Headers/simulation/director.h
        Headers/simulation/simulation.h
        Headers/simulation/simulation2dap.h
        Headers/simulation/simulation3dcs.h
        Headers/ui/menu.h
        Headers/utils/datasaver.h
        Headers/utils/leapfrog.h
        Headers/utils/progressbar.h
        Headers/utils/vector_operations.h
)

# Main executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} inih matplot)


add_subdirectory(Google_tests)
