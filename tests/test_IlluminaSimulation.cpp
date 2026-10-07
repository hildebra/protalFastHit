// The Illumina paired-end reads simulate_metagenomes makes in process (src/RandomForest/IlluminaSimulator.h, on
// ReadPipeline.h): the instrument profiles against their sources, the pairs from their fragments, names, a host's
// pairs, the same files for any number of threads, first reads alone, and samples streamed into named pipes.

#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <future>
#include <map>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <sys/stat.h>

#include "IO/ThreadedGzStream.h"
#include "RandomForest/IlluminaSimulator.h"
#include "TestUtil.h"

using namespace protal::sim;
using protal::test::ScratchDir;

namespace {

struct Fq {
    std::string name, seq, qual;
};

std::vector<Fq> ReadFastq(std::string const& path) {
    protal::ThreadedGzIstream in(path.c_str());
    EXPECT_TRUE(in.rdbuf()->is_open()) << path;
    std::vector<Fq> reads;
    std::string name, seq, plus, qual;
    while (std::getline(in, name) && std::getline(in, seq) && std::getline(in, plus) && std::getline(in, qual)) {
        reads.push_back({name, seq, qual});
    }
    return reads;
}

std::string Revcomp(std::string s) {
    ReverseComplement(s);
    return s;
}

int Mismatches(std::string const& a, std::string const& b) {
    int m = 0;
    for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i) m += a[i] != b[i];
    return m + static_cast<int>(std::max(a.size(), b.size()) - std::min(a.size(), b.size()));
}

IlluminaSetup MakeSetup(std::string const& profile, int length, double fragment = 350) {
    IlluminaSetup setup;
    setup.profile = IlluminaProfile::Named(profile);
    setup.read_length = length;
    setup.fragment_mean = fragment;
    setup.fragment_sd = 50;
    return setup;
}

}  // namespace

TEST(IlluminaSimulation, ProfilesMatchTheirSources) {
    auto const names = IlluminaProfile::Names();
    for (auto const* name : {"HS20", "HS25", "HSXt", "NovaSeq", "MSv3"}) {
        EXPECT_NE(std::find(names.begin(), names.end(), name), names.end()) << name;
    }
    EXPECT_THROW(IlluminaProfile::Named("HiSeq"), std::invalid_argument);

    // NovaSeq 6000: qualities written in 4 bins, ~94% at Q37 and a few % at 23 and 12 (Illumina's RTA3 note); errors
    // R1 ~0.14-0.35% and R2 ~0.2-0.6% (PhiX), R2 the worse.
    auto const nova = ReportProfile(MakeSetup("NovaSeq", 150), 4000);
    for (int r = 0; r < 2; ++r) {
        for (int q = 0; q < 94; ++q) {
            if (nova.quality_share[r][q] > 0) EXPECT_TRUE(q == 2 || q == 12 || q == 23 || q == 37) << q;
        }
    }
    double const q37 = (nova.quality_share[0][37] + nova.quality_share[1][37]) / 2;
    EXPECT_NEAR(q37, 0.94, 0.03);
    EXPECT_GT(nova.substitutions[0], 0.0010);
    EXPECT_LT(nova.substitutions[0], 0.0040);
    EXPECT_GT(nova.substitutions[1], nova.substitutions[0]);
    EXPECT_LT(nova.substitutions[1], 0.0070);
    EXPECT_LT(nova.insertions[0] + nova.deletions[0], 5e-5);

    // HiSeq 2500: the 8-level bins only.
    auto const hs25 = ReportProfile(MakeSetup("HS25", 125), 2000);
    std::set<int> const bins = {2, 6, 15, 22, 27, 33, 37, 40};
    for (int q = 0; q < 94; ++q) {
        if (hs25.quality_share[0][q] > 0) EXPECT_TRUE(bins.count(q)) << q;
    }

    // MiSeq v3 2x300: read 2 falls to ~Q17 at its end, read 1 to ~Q23 (InSilicoSeq's MiSeq curve); >70% of the bases
    // at Q30 or more (the specification); insertions ~4e-5 per base (Schirmer 2015).
    auto const miseq = ReportProfile(MakeSetup("MSv3", 300, 550), 2000);
    EXPECT_NEAR(miseq.mean_by_cycle[0].back(), 23.0, 1.5);
    EXPECT_NEAR(miseq.mean_by_cycle[1].back(), 17.0, 1.5);
    EXPECT_GT((miseq.q30[0] + miseq.q30[1]) / 2, 0.68);
    EXPECT_NEAR(miseq.insertions[0], 4e-5, 2e-5);
    EXPECT_GT(miseq.substitutions[1], 2 * miseq.substitutions[0]);

    // The calibrated curves follow the published ones; --mean_quality shifts them (and the errors with them).
    auto const hsxt = ReportProfile(MakeSetup("HSXt", 150), 2000);
    IlluminaModel const model(MakeSetup("HSXt", 150));
    for (int c : {0, 74, 149}) EXPECT_NEAR(hsxt.mean_by_cycle[0][c], model.TargetMean(0)[c], 1.0) << c;
    auto shifted = MakeSetup("HSXt", 150);
    shifted.mean_quality = 30;
    auto const q30 = ReportProfile(shifted, 2000);
    for (int r = 0; r < 2; ++r) {
        double mean = 0;
        for (double q : q30.mean_by_cycle[r]) mean += q;
        EXPECT_NEAR(mean / 150, 30.0, 1.0) << r;
        EXPECT_GT(q30.substitutions[r], hsxt.substitutions[r]);
    }
}

TEST(IlluminaSimulation, PairsOfTheirFragments) {
    std::mt19937 gen(5);
    std::string const genome = protal::test::RandomSequence(20000, gen);
    IlluminaModel const model(MakeSetup("NovaSeq", 150));
    LongRng rng(3);
    MadeRead r1, r2;
    int close = 0, both_strands = 0;
    for (int i = 0; i < 400; ++i) {
        std::uint32_t const length = model.FragmentLength(rng, genome.size());
        ASSERT_GE(length, 150u);
        std::string const fragment = genome.substr(rng.Below(genome.size() - length + 1), length);
        model.Pair(fragment, 0.0, rng, r1, r2);
        ASSERT_EQ(r1.seq.size(), 150u);
        ASSERT_EQ(r2.seq.size(), 150u);
        ASSERT_EQ(r1.qual.size(), 150u);
        // read 1: the fragment's start; read 2: its reverse complement's start (a few errors, rarely an indel)
        bool const ok1 = Mismatches(r1.seq, fragment.substr(0, 150)) <= 6;
        bool const ok2 = Mismatches(r2.seq, Revcomp(fragment).substr(0, 150)) <= 6;
        close += ok1 && ok2;
        both_strands += r1.seq != r2.seq;
    }
    EXPECT_GE(close, 380);
    EXPECT_EQ(both_strands, 400);
    // fragment lengths: normal (350, 50), at least a read long
    double sum = 0;
    for (int i = 0; i < 4000; ++i) sum += model.FragmentLength(rng, 100000);
    EXPECT_NEAR(sum / 4000, 350, 5);
    EXPECT_EQ(model.FragmentLength(rng, 100), 150u);  // contigs shorter than a read: a read long
}

TEST(IlluminaSimulation, SamplesNamesHostThreadsAndPipes) {
    ScratchDir dir("illumina");
    std::mt19937 gen(9);
    std::string const a1 = protal::test::RandomSequence(6000, gen), a2 = protal::test::RandomSequence(120, gen);
    std::string const b = protal::test::RandomSequence(30000, gen);
    auto const ga = dir.Write("GA.fna", ">contigA1 x\n" + a1 + "\n>contigA2\n" + a2 + "\n");
    auto const gb = dir.Write("GB.fna", ">contigB\n" + b + "\n");
    std::string const h1 = protal::test::RandomSequence(15000, gen);
    std::filesystem::create_directories(dir.path / "host");
    dir.Write("host/host.seq", h1);
    dir.Write("host/host.json", "{\"source\": [\"/x\", 1, 2], \"contigs\": [[\"chr1\", 0, 15000]], \"bases\": 15000}");

    auto samples_in = [&](std::filesystem::path const& out, std::string const& suffix) {
        std::vector<PairedSample> samples(3);
        samples[0] = {"s1", out / ("s1_R1" + suffix), out / ("s1_R2" + suffix), {{"GA", ga, 1000, 11}, {"GB", gb, 2500, 12}}};
        samples[1] = {"s2", out / ("s2_R1" + suffix), out / ("s2_R2" + suffix), {{"GB", gb, 700, 13}}};
        samples[1].host_pairs = 300;
        samples[1].host_seed = 77;
        samples[2] = {"s3", out / ("s3_R1" + suffix), "", {{"GA", ga, 1000, 11}, {"GB", gb, 2500, 12}}};  // first reads
        for (auto& s : samples) s.run_seed = 5;
        return samples;
    };
    PairedOptions options;
    options.setup = MakeSetup("NovaSeq", 150);
    options.host = dir.path / "host";
    options.threads = 1;
    auto const one = samples_in(dir.path / "one", ".fq.zst");
    auto const totals = SimulatePairs(one, options);
    options.threads = 4;
    auto const four = samples_in(dir.path / "four", ".fq.zst");
    SimulatePairs(four, options);
    auto const gz = samples_in(dir.path / "gz", ".fq.gz");
    SimulatePairs(gz, options);
    ASSERT_EQ(totals.size(), 3u);
    EXPECT_EQ(totals[0].reads, 3500u);
    EXPECT_EQ(totals[1].reads, 1000u);
    for (std::size_t s = 0; s < 3; ++s) {
        EXPECT_EQ(protal::test::Slurp(one[s].r1), protal::test::Slurp(four[s].r1)) << "same files on 1 and 4 threads";
        if (!one[s].r2.empty()) EXPECT_EQ(protal::test::Slurp(one[s].r2), protal::test::Slurp(four[s].r2));
        auto const zst = ReadFastq(one[s].r1.string()), bgzf = ReadFastq(gz[s].r1.string());
        ASSERT_EQ(zst.size(), bgzf.size());
        for (std::size_t i = 0; i < zst.size(); ++i) EXPECT_EQ(zst[i].seq, bgzf[i].seq);
    }
    // first reads alone: the same as the R1 of a run with both
    EXPECT_FALSE(std::filesystem::exists(dir.path / "one" / "s3_R2.fq.zst"));
    auto const s1r1 = ReadFastq(one[0].r1.string()), s3r1 = ReadFastq(one[2].r1.string());
    ASSERT_EQ(s1r1.size(), s3r1.size());
    for (std::size_t i = 0; i < s1r1.size(); ++i) ASSERT_EQ(s1r1[i].seq, s3r1[i].seq);

    // names <contig>-<n>/1 and /2, n = 1, 2, ... over the sample; reads from their contig (contigA2 is shorter than a
    // read: none); the host's h_<n>
    auto const r1 = ReadFastq(one[0].r1.string()), r2 = ReadFastq(one[0].r2.string());
    ASSERT_EQ(r1.size(), 3500u);
    ASSERT_EQ(r2.size(), 3500u);
    std::map<std::string, int> per_contig;
    int located = 0;
    for (std::size_t i = 0; i < r1.size(); ++i) {
        auto const& name = r1[i].name;
        auto const dash = name.rfind('-');
        ASSERT_NE(dash, std::string::npos);
        EXPECT_EQ(name.substr(dash + 1), std::to_string(i + 1) + "/1");
        EXPECT_EQ(r2[i].name, name.substr(0, name.size() - 1) + "2");
        std::string const contig = name.substr(1, dash - 1);
        ++per_contig[contig];
        std::string const& seq = contig == "contigA1" ? a1 : b;
        std::string const head = r1[i].seq.substr(0, 25);
        located += seq.find(head) != std::string::npos || seq.find(Revcomp(head)) != std::string::npos;
    }
    EXPECT_EQ(per_contig["contigA1"], 1000);
    EXPECT_EQ(per_contig["contigB"], 2500);
    EXPECT_EQ(per_contig.count("contigA2"), 0u);
    EXPECT_GE(located, 3300);
    auto const host = ReadFastq(one[1].r1.string());
    ASSERT_EQ(host.size(), 1000u);
    int host_reads = 0, in_host = 0;
    for (auto const& read : host) {
        if (read.name.rfind("@h_", 0) != 0) continue;
        ++host_reads;
        std::string const head = read.seq.substr(0, 25);
        in_host += h1.find(head) != std::string::npos || h1.find(Revcomp(head)) != std::string::npos;
    }
    EXPECT_EQ(host_reads, 300);
    EXPECT_GE(in_host, 280);
    EXPECT_EQ(host.back().name, "@h_1000/1");  // after the community's 700

    // Streamed: the outputs named pipes, read sample by sample, R1 and R2 in step (32 records of each in turn, as
    // protal reads them): the same reads as the files.
    std::filesystem::create_directories(dir.path / "pipes");
    auto piped = samples_in(dir.path / "pipes", ".fq.zst");
    piped.resize(2);
    for (auto const& s : piped) {
        ASSERT_EQ(::mkfifo(s.r1.c_str(), 0600), 0);
        ASSERT_EQ(::mkfifo(s.r2.c_str(), 0600), 0);
    }
    auto reader = std::async(std::launch::async, [&] {
        std::vector<std::vector<std::string>> got;
        for (auto const& s : piped) {
            protal::ThreadedGzIstream in1(s.r1.c_str());
            protal::ThreadedGzIstream in2(s.r2.c_str());
            std::vector<std::string> seqs;
            std::string line[4];
            bool more = true;
            while (more) {
                for (int k = 0; k < 32 && more; ++k) {
                    for (auto& l : line) more = more && static_cast<bool>(std::getline(in1, l));
                    if (more) seqs.push_back(line[1]);
                }
                for (int k = 0; k < 32; ++k) {
                    bool ok = true;
                    for (auto& l : line) ok = ok && static_cast<bool>(std::getline(in2, l));
                    if (!ok) break;
                    seqs.push_back(line[1]);
                }
            }
            got.push_back(seqs);
        }
        return got;
    });
    options.threads = 3;
    SimulatePairs(piped, options);
    ASSERT_EQ(reader.wait_for(std::chrono::seconds(60)), std::future_status::ready);
    auto const got = reader.get();
    ASSERT_EQ(got.size(), 2u);
    for (std::size_t s = 0; s < 2; ++s) {
        auto const f1 = ReadFastq(one[s].r1.string()), f2 = ReadFastq(one[s].r2.string());
        EXPECT_EQ(got[s].size(), f1.size() + f2.size()) << s;
        std::multiset<std::string> expected, seen(got[s].begin(), got[s].end());
        for (auto const& r : f1) expected.insert(r.seq);
        for (auto const& r : f2) expected.insert(r.seq);
        EXPECT_TRUE(expected == seen) << s;
    }
}
