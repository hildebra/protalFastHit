find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(protal_tests
        test_SNPCalling.cpp
        test_SamRoundTrip.cpp
        test_Operations.cpp
        test_Index.cpp
        test_InputValidation.cpp
        test_Parsing.cpp
        test_StrainOutput.cpp
        test_Zstd.cpp
        test_IndexCodec.cpp
        test_Database.cpp
        test_LongReads.cpp)
target_compile_options(protal_tests PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_libraries(protal_tests PRIVATE protal_lib dl z isa_baseline GTest::gtest_main)

gtest_discover_tests(protal_tests)
