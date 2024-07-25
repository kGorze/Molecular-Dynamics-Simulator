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

add_subdirectory(Third_party/matplotplusplus)



# inih library
add_library(inih STATIC third_party/inih/ini.c)
# target_include_directories(inih PUBLIC third_party/inih)
# target_compile_features(inih PUBLIC cxx_std_23)

# Source files
set(SOURCES
        main.cpp
        simulation.cpp
        builder.cpp
        director.cpp
        abstract_factory.cpp
        menu.cpp
        concretecommands.cpp
        command.cpp
        Headers/vector_operations.h
        Headers/types.h
        Headers/leapfrog.h
        Headers/boundary_conditions.h
        Headers/progressbar.h
        Headers/simulation2dap.h
        Headers/simulation3dcs.h
        progressbar.cpp
        errorcodes.cpp
        leapfrog.cpp
        simulation2dap.cpp
        simulation3dcs.cpp
        Headers/datasaver.h
        datasaver.cpp
        Headers/constants.h
)

# Main executable
add_executable(${PROJECT_NAME} ${SOURCES})
# add_executable(${PROJECT_NAME} ${SOURCES} ${ZMQ})
# target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_23)

# Link libraries
target_link_libraries(${PROJECT_NAME} inih matplot)
# target_link_libraries(${PROJECT_NAME} stdc++ zmq)
