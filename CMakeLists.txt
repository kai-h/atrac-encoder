cmake_minimum_required(VERSION 3.16)
project(atrac3 VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
# Self-contained executables: static C runtime on Windows
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

add_library(atrac3_lib STATIC
  src/tables.c
  src/lp2_analysis.c src/lp2_alloc.c src/lp2_bitstream.c
  src/lp4_analysis.c src/lp4_gain.c src/lp4_stereo.c src/lp4_alloc.c src/lp4_encoder.c
  src/encoder.c src/decoder.c src/files.c)
set_target_properties(atrac3_lib PROPERTIES OUTPUT_NAME atrac3)
target_include_directories(atrac3_lib PUBLIC src)

# Byte-exact output needs every float/double operation to round exactly as written:
# no fused multiply-add, no excess precision (x87), no fast-math reassociation.
if(MSVC)
  target_compile_options(atrac3_lib PRIVATE /fp:precise)
else()
  target_compile_options(atrac3_lib PRIVATE -ffp-contract=off -fno-fast-math)
  if(NOT EMSCRIPTEN AND (CMAKE_SYSTEM_PROCESSOR MATCHES "^(i.86|x86)$" OR
     (CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$" AND CMAKE_SIZEOF_VOID_P EQUAL 4)))
    target_compile_options(atrac3_lib PRIVATE -msse2 -mfpmath=sse)
  endif()
endif()
if(UNIX AND NOT APPLE AND NOT EMSCRIPTEN)
  target_link_libraries(atrac3_lib PUBLIC m)
endif()

if(EMSCRIPTEN)
  # Browser module: ES6 factory exporting the whole-file API and the streaming encoder (see web/)
  add_executable(atrac3_wasm web/glue.c)
  set_target_properties(atrac3_wasm PROPERTIES OUTPUT_NAME atrac3-core SUFFIX ".mjs")
  target_link_libraries(atrac3_wasm atrac3_lib)
  target_link_options(atrac3_wasm PRIVATE
    -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createAtrac3Module
    -sENVIRONMENT=web,worker,node -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=1048576
    "-sEXPORTED_FUNCTIONS=['_at3w_encode','_at3w_decode','_at3w_result','_at3w_result_len','_at3_encoder_create','_at3_encoder_destroy','_at3_encoder_frame_bytes','_at3_encoder_max_output','_at3_encoder_write','_at3_encoder_finish','_malloc','_free']"
    "-sEXPORTED_RUNTIME_METHODS=['HEAPU8','UTF8ToString']")
  configure_file(web/atrac3.mjs atrac3.mjs COPYONLY)
  configure_file(web/index.html index.html COPYONLY)
else()
  add_executable(atrac3 cli/atrac3.c)
  target_link_libraries(atrac3 atrac3_lib)
  target_compile_definitions(atrac3 PRIVATE AT3_VERSION="${PROJECT_VERSION}")
  option(ATRAC3_STATIC "Link the command-line tool statically (Linux release builds)" OFF)
  if(ATRAC3_STATIC)
    target_link_options(atrac3 PRIVATE -static)
  endif()
  install(TARGETS atrac3 RUNTIME DESTINATION bin)
endif()
