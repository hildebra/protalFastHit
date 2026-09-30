#include "probe.cpp"

// Q3: reads with a deletion over a position are counted as its reference allele.
TEST(Probe2, DeletionReadsBecomeReferenceSupportAtASnp) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    char const alt = Alt(reference[21]);
    auto snp = reference; snp[21] = alt;
    std::string del = reference.substr(0, 20) + reference.substr(23);  // 3 bp deletion 20-22
    for (double af : { 0.0, 0.2 }) {
        StrainLevelContainer s(ref.gene());
        for (int i = 0; i < 8; i++) ASSERT_TRUE(s.AddSam(MakeSam(snp, "21M1X28M", 1, i % 2 ? 0x10 : 0), i, true));
        for (int i = 0; i < 2; i++) ASSERT_TRUE(s.AddSam(MakeSam(del, "20M3D27M", 1, i % 2 ? 0x10 : 0), 10 + i, true));
        s.PostProcess(2, 2, af, 15, 90, true);
        std::cout << "af " << af << " bin21: " << BinString(s, 21) << '\n';
        auto r = RunMSA({ &s }, reference, { af });
        Print("8 reads SNP at 21, 2 reads with a deletion over 20-22, af " + std::to_string(af), r);
        EXPECT_EQ(r.rows[0][21], alt) << "no read shows the reference base at 21";
    }
}
