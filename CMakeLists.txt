cmake_minimum_required(VERSION 3.20)

project(predicates VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/CPM.cmake)

CPMAddPackage(
  NAME SIMDE
  URL https://github.com/simd-everywhere/simde/archive/refs/tags/v0.8.4-rc2.zip
)

add_library(simde INTERFACE)
target_include_directories(simde INTERFACE "${SIMDE_SOURCE_DIR}/simde")
target_compile_definitions(simde INTERFACE SIMDE_ENABLE_NATIVE_ALIASES)
if(APPLE)
  target_compile_options(simde INTERFACE -march=armv8-a+simd)
endif()

add_library(predicates STATIC
  src/predicates.cpp
)

target_include_directories(predicates PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(predicates PUBLIC simde)
target_compile_features(predicates PRIVATE cxx_std_23)

add_executable(predicates_cli
  src/main.cpp
)

target_link_libraries(predicates_cli PRIVATE
  predicates
)

include(CTest)
if(BUILD_TESTING)
  add_executable(predicates_tests
    tests/predicates_test.cpp
  )
  target_link_libraries(predicates_tests PRIVATE predicates)
  add_test(NAME predicates_tests COMMAND predicates_tests)
endif()
