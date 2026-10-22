
cmake_minimum_required(VERSION 3.20)
project(GAMES201_MPM88 LANGUAGES CXX)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EIGEN_ROOT "D:/games101/eigen-3.4.0"
    CACHE PATH "Directory containing Eigen/Dense")

set(OpenCV_DIR "D:/games101/opencv/build"
    CACHE PATH "Directory containing OpenCVConfig.cmake")
if(NOT EXISTS "${EIGEN_ROOT}/Eigen/Dense")
    message(FATAL_ERROR
        "Eigen3 not found! Please check EIGEN_ROOT.")
endif()
find_package(OpenCV REQUIRED COMPONENTS
    core
    highgui
    imgproc
    imgcodecs
)
add_executable(mpm88
    mpm88.cpp
)

target_include_directories(mpm88 SYSTEM PRIVATE
    "${EIGEN_ROOT}"
    ${OpenCV_INCLUDE_DIRS}
)
target_link_libraries(mpm88 PRIVATE
    ${OpenCV_LIBS}
)


if(MSVC)
    target_compile_options(mpm88 PRIVATE
        /utf-8
        /W4
        /external:W0
    )
    target_compile_options(mpm88 PRIVATE
        $<$<CONFIG:Release>:/O2>
    )
endif()
foreach(lib IN LISTS OpenCV_LIBS)
    if(TARGET ${lib})
        get_target_property(lib_type ${lib} TYPE)

        if(lib_type STREQUAL "SHARED_LIBRARY")
            add_custom_command(TARGET mpm88 POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "$<TARGET_FILE:${lib}>"
                    "$<TARGET_FILE_DIR:mpm88>"
                VERBATIM
            )
        endif()
    endif()
endforeach()

set_property(DIRECTORY PROPERTY
    VS_STARTUP_PROJECT mpm88
)
