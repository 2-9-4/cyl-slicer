cmake_minimum_required(VERSION 3.14)
project(gcode_viewer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)

# GLFW
find_package(glfw3 QUIET)
if(NOT glfw3_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLFW REQUIRED glfw3)
endif()

# GLEW
find_package(GLEW QUIET)
if(NOT GLEW_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLEW REQUIRED glew)
endif()

# ImGui via FetchContent
include(FetchContent)
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.6
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(imgui)
if(NOT imgui_POPULATED)
    FetchContent_Populate(imgui)
    add_library(imgui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
    )
    target_include_directories(imgui PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
    )
    if(glfw3_FOUND)
        target_link_libraries(imgui PUBLIC glfw OpenGL::GL)
    else()
        target_include_directories(imgui PUBLIC ${GLFW_INCLUDE_DIRS})
        target_link_libraries(imgui PUBLIC ${GLFW_LIBRARIES} OpenGL::GL)
    endif()
endif()

add_executable(gcode_viewer gcode_viewer.cpp)
target_compile_options(gcode_viewer PRIVATE -O2 -Wall)
target_link_libraries(gcode_viewer imgui)

if(glfw3_FOUND AND GLEW_FOUND)
    target_link_libraries(gcode_viewer OpenGL::GL glfw GLEW::GLEW)
elseif(glfw3_FOUND)
    target_include_directories(gcode_viewer PRIVATE ${GLEW_INCLUDE_DIRS})
    target_link_libraries(gcode_viewer OpenGL::GL glfw ${GLEW_LIBRARIES})
else()
    target_include_directories(gcode_viewer PRIVATE ${GLFW_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
    target_link_libraries(gcode_viewer OpenGL::GL ${GLFW_LIBRARIES} ${GLEW_LIBRARIES})
endif()
