include(FetchContent)

FetchContent_Declare(
  scran_tests
  GIT_REPOSITORY https://github.com/libscran/scran_tests
  GIT_TAG master
)

FetchContent_MakeAvailable(scran_tests)

add_executable(
    libtest 
    src/auc.cpp
    src/cohens_d.cpp
    src/simple_diff.cpp
    src/summarize_comparisons.cpp
    src/summarize_effects.cpp
    src/score_markers_pairwise.cpp
    src/score_markers_summary.cpp
    src/create_combinations.cpp
    src/average_group_stats.cpp
)

target_link_libraries(
    libtest
    scran_markers
    scran_tests
)

target_compile_options(libtest PRIVATE -Wall -Werror -Wpedantic -Wextra)

option(CODE_COVERAGE "Enable coverage testing" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

include(GoogleTest)
gtest_discover_tests(libtest)
