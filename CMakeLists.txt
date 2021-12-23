cmake_minimum_required(VERSION 3.16)

PROJECT (Lab1)

add_definitions(-DCMAKE_BUILD_TYPE=Debug)

#头文件目录
# include_directories("D:\\2021Fall\\ComputerGraphics\\openGL4.6\\Include")
include_directories(./include)
include_directories("D:\\2021Fall\\ComputerGraphics\\openGL4.6\\Include")
#库文件目录
link_directories("D:\\2021Fall\\ComputerGraphics\\openGL4.6\\Libs")

#项目文件夹下源文件有main.cpp
# aux_source_directory(./src DIR_ALL_SRC)				   # src目录下所有文件命名为DIR_ALL_SRC
# Add project files
set(SRC_FILES ./src/main.cpp
              ./src/CageGenerator.cpp
              ./src/Geometry.cpp
              ./src/glad.c
              ./src/glRender.cpp)

add_compile_options(-Wa,-mbig-obj)#  “File too big/too many sections” error. => set the assembler’s flags by using -Wa,<assembler_flags_comma_separated> in compiler
add_compile_options(-Os)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,1000000000")

add_executable(Lab1 ${SRC_FILES})

target_link_libraries(Lab1 glfw3.dll)