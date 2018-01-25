cmake_minimum_required(VERSION 2.6)
project(GcodePovray)

aux_source_directory(. src)
add_executable(GcodePovray ${src})