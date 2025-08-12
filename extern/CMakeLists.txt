include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # ^3.0.0
)

FetchContent_Declare(
  tatami_stats
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_stats
  GIT_TAG master # ^1.0.0
)

FetchContent_Declare(
  scran_blocks
  GIT_REPOSITORY https://github.com/libscran/scran_blocks
  GIT_TAG master # ^0.1.0
)

FetchContent_Declare(
  sanisizer
  GIT_REPOSITORY https://github.com/LTLA/sanisizer
  GIT_TAG master # ^0.1.1
)

FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(tatami_stats)
FetchContent_MakeAvailable(scran_blocks)
FetchContent_MakeAvailable(sanisizer)
