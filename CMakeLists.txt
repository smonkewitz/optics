cmake_minimum_required(VERSION 3.24)

project(
  optics
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

set(ABSL_PROPAGATE_CXX_STD ON)

add_compile_options(-Wall -Wextra)

include(FetchContent)

# Disable installation of embedded googletest
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

# googletest v1.14.0
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/f8d7d77c06936315286eb55f8de22cd23c188571.zip
)
FetchContent_MakeAvailable(googletest)

# Abseil 20240116.0
FetchContent_Declare(
  absl
  URL https://github.com/abseil/abseil-cpp/archive/4a2c63365eff8823a5221db86ef490e828306f9d.zip
)
FetchContent_MakeAvailable(absl)

# fmtlib 10.2.1
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/e69e5f977d458f2650bb346dadf2ad30c5320281.zip
)
FetchContent_MakeAvailable(fmt)

# fast_float 6.0.0
FetchContent_Declare(
  fast_float
  URL https://github.com/fastfloat/fast_float/archive/f3ff46fd40d28b5eb272727812e5bd8dc5a6d876.zip
)
FetchContent_MakeAvailable(fast_float)

enable_testing()
add_subdirectory(src)
