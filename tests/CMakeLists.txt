find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(protal_tests
        test_SNPCalling.cpp
        test_SamRoundTrip.cpp
        test_RunStatusAndBgzf.cpp
        test_Index.cpp
        test_PackedIndex.cpp
        test_InputValidation.cpp
        test_Parsing.cpp
        test_StrainOutput.cpp
        test_SequenceRange.cpp
        test_Zstd.cpp
        test_IndexCodec.cpp
        test_Database.cpp
        test_LongReads.cpp
        test_ReaderAndTimers.cpp
        test_WFA2Wrapper.cpp
        test_Syncmers.cpp
        test_AnchoredAlignment.cpp
        test_SamFile.cpp
        test_CommunityDesign.cpp
        test_LongReadSimulation.cpp
        test_IlluminaSimulation.cpp
        test_GenomeStore.cpp
        test_PackedSequence.cpp
        test_ReverseComplement.cpp
        test_GeneWindows.cpp
        test_PairJoin.cpp
        test_GeneNeighbours.cpp
        test_Haplotypes.cpp
        test_ProfileThreads.cpp
        test_SampleContext.cpp
        test_FalsePositiveFeatures.cpp
        test_ComplexCommunityFeatures.cpp
        test_ReadEvidence.cpp
        test_AlignmentScreen.cpp
        test_GeneTableFile.cpp
        test_ExactSum.cpp
        test_FlexScan.cpp
        test_SharedSeeds.cpp
        test_GoldenModelRules.cpp
        test_GeneConservation.cpp
        test_CongenerGaps.cpp
        test_StrainAlleles.cpp
        test_ModelFeatures.cpp
        test_AncestrySites.cpp
        test_Composition.cpp)
target_compile_options(protal_tests PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
# tests/data (golden_model_rules.tsv), wherever the tests run from; on that file alone, so the others' compile commands
# (and their compiler cache entries) stay as they were.
set_source_files_properties(test_GoldenModelRules.cpp PROPERTIES COMPILE_DEFINITIONS PROTAL_TESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(protal_tests PRIVATE protal_lib dl isa_baseline GTest::gtest_main)

gtest_discover_tests(protal_tests)
