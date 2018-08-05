cmake_minimum_required(VERSION 2.8)
project(tinymp3)
add_executable(tinymp3 main.c shine_mp3.c)
target_link_libraries(tinymp3 -lm)