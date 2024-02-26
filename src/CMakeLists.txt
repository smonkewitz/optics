add_library(optics-lib)

target_sources(
  optics-lib
  PRIVATE
    InputFile.cc
    LonLat.cc
    Optics.cc
    SeedList.cc
    Tree.cc
)

target_link_libraries(
  optics-lib
  PUBLIC
    absl::check
    absl::log
    fast_float
    fmt::fmt
)

# Unit tests

add_executable(optics-test)

target_sources(
  optics-test
  PRIVATE
    TreeTest.cc
    SeedListTest.cc
)

target_link_libraries(
  optics-test
  PUBLIC
    optics-lib
    GTest::gmock_main
)

include(GoogleTest)

gtest_discover_tests(optics-test)
