include(FetchContent)

FetchContent_Declare(
  scran_tests
  GIT_REPOSITORY https://github.com/libscran/scran_tests
  GIT_TAG master
)

FetchContent_MakeAvailable(scran_tests)

include(GoogleTest)

option(CODE_COVERAGE "Enable coverage testing" OFF)
set(DO_CODE_COVERAGE OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(DO_CODE_COVERAGE ON)
endif()

macro(decorate_test target)
    target_link_libraries(
        ${target}
        scran_markers
        scran_tests
    )

    target_compile_options(${target} PRIVATE -Wall -Werror -Wpedantic -Wextra)

    if(DO_CODE_COVERAGE)
        target_compile_options(${target} PRIVATE -O0 -g --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endif()

    gtest_discover_tests(${target})
endmacro()

set(LIST_FILES
    src/auc.cpp
    src/cohens_d.cpp
    src/simple_diff.cpp
    src/summarize_comparisons.cpp
    src/summarize_effects.cpp
    src/score_markers_pairwise.cpp
    #src/score_markers_summary.cpp
    #src/score_markers_best.cpp
    src/create_combinations.cpp
    src/average_group_stats.cpp
)

add_executable(libtest "${LIST_FILES}")
decorate_test(libtest)

add_executable(dirtytest "${LIST_FILES}")
decorate_test(dirtytest)
target_compile_definitions(dirtytest PRIVATE "SCRAN_MARKERS_TEST_INIT=scran_tests::initial_value()")
