// The long and Ultima reads simulate_metagenomes makes in process (src/RandomForest/LongReadSimulator.h): the read
// models (hifi_reads.py's HiFi and flow models, pbsim3's qshmm), the templates drawn from genomes and a host, the
// reads' names and files, and that they do not depend on the threads.

#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "IO/ThreadedGzStream.h"
#include "RandomForest/LongReadSimulator.h"
#include "TestUtil.h"

using namespace protal::sim;
using protal::test::ScratchDir;

namespace {

std::string Revcomp(std::string s) {
    ReverseComplement(s);
    return s;
}

// The reads of a compressed FASTQ (zstd or BGZF): {name, sequence, qualities}.
struct FastqRead {
    std::string name, seq, qual;
};

std::vector<FastqRead> ReadFastq(std::string const& path) {
    protal::ThreadedGzIstream in(path.c_str());
    EXPECT_TRUE(in.rdbuf()->is_open()) << path;
    std::vector<FastqRead> reads;
    std::string name, seq, plus, qual;
    while (std::getline(in, name) && std::getline(in, seq) && std::getline(in, plus) && std::getline(in, qual)) {
        reads.push_back({name, seq, qual});
    }
    EXPECT_FALSE(in.rdbuf()->read_failed()) << path;
    return reads;
}

double Phred(char q) { return q - 33; }

}  // namespace

TEST(LongReadSimulation, RandomStreamsAndSetups) {
    LongRng a(5), b(5), c(6);
    std::vector<std::uint64_t> sa, sb, sc;
    for (int i = 0; i < 4; ++i) {
        sa.push_back(a());
        sb.push_back(b());
        sc.push_back(c());
    }
    EXPECT_EQ(sa, sb);
    EXPECT_NE(sa, sc);
    std::vector<int> counts(7, 0);
    double sum = 0, squares = 0;
    for (int i = 0; i < 70000; ++i) {
        ++counts[a.Below(7)];
        double const x = a.Normal();
        sum += x;
        squares += x * x;
    }
    for (int n : counts) EXPECT_NEAR(n, 10000, 400);
    EXPECT_NEAR(sum / 70000, 0.0, 0.02);
    EXPECT_NEAR(squares / 70000, 1.0, 0.03);
    EXPECT_NE(MixSeed(1, 2), MixSeed(2, 1));

    auto hifi = LongReadSetup::Parse("hifi:15000:3000:3");
    EXPECT_EQ(hifi.method, LongReadSetup::Method::Hifi);
    EXPECT_EQ(hifi.length_mean, 15000);
    EXPECT_EQ(hifi.q_sd, 3);
    auto ultima = LongReadSetup::Parse("ultima:300:40:25:2");
    EXPECT_EQ(ultima.method, LongReadSetup::Method::Ultima);
    EXPECT_EQ(ultima.q_mean, 25);
    auto ont = LongReadSetup::Parse("qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36");
    EXPECT_EQ(ont.method, LongReadSetup::Method::Qshmm);
    EXPECT_EQ(ont.model, "QSHMM-ONT-HQ");
    EXPECT_EQ(ont.accuracy, 0.97);
    EXPECT_EQ(ont.sub_ratio, 39);
    EXPECT_EQ(ont.del_ratio, 36);
    EXPECT_EQ(LongReadSetup::Parse("qshmm:M:8000:6000:0.9").ins_ratio, 55);  // pbsim3's default ratio 6:55:39
    for (auto const* bad : {"errhmm:ERRHMM-SEQUEL:15000:3000:0.99", "hifi:15000:3000", "qshmm:M:8000:6000:1.5",
                            "qshmm:M:8000:6000:0.9:1/2", "nanopore:1:2:3", "hifi:x:3000:3", "hifi:0:3000:3"}) {
        EXPECT_THROW(LongReadSetup::Parse(bad), std::invalid_argument) << bad;
    }

    // A read's length: gamma of the mean and SD within [100, 1 Mb]; the mean if the SD is 0.
    double length_sum = 0;
    for (int i = 0; i < 20000; ++i) {
        auto const length = DrawReadLength(a, 8000, 6000);
        ASSERT_GE(length, 100u);
        ASSERT_LE(length, 1000000u);
        length_sum += length;
    }
    EXPECT_NEAR(length_sum / 20000, 8000, 150);
    EXPECT_EQ(DrawReadLength(a, 1000, 0), 1000u);
    EXPECT_EQ(DrawReadLength(a, 50, 0), 100u);

    EXPECT_EQ(Revcomp("AACGTNRx"), "xRNACGTT");
}

TEST(LongReadSimulation, HifiAndFlowModels) {
    // hifi_reads.py's read quality by length.
    std::vector<double> q;
    for (double length : {1000.0, 5000.0, 15000.0, 25000.0, 37500.0, 50000.0, 75000.0}) q.push_back(hifi::LengthQuality(length));
    EXPECT_EQ(q, (std::vector<double>{50, 50, 40, 30, 25, 20, 10}));

    // HiFi reads: their qualities say how many errors they have (the expected errors of the qualities against the
    // errors made), mostly in homopolymers, and the read's mean error rate is about its length's quality.
    std::mt19937 gen(3);
    auto const setup = LongReadSetup::Parse("hifi:15000:0:3");
    LongRng rng(11);
    MadeRead read;
    std::uint64_t errors = 0;
    double expected = 0, length_sum = 0;
    for (int i = 0; i < 300; ++i) {
        std::string const templ = protal::test::RandomSequence(25000, gen);
        hifi::Mutate(templ, rng, setup, read, &errors);
        ASSERT_EQ(read.seq.size(), read.qual.size());
        length_sum += read.seq.size();
        for (char c : read.qual) {
            ASSERT_GE(Phred(c), 1);
            ASSERT_LE(Phred(c), 93);
            expected += std::pow(10.0, -Phred(c) / 10.0);
        }
    }
    EXPECT_NEAR(errors / expected, 1.0, 0.1) << errors << " errors, " << expected << " expected by the qualities";
    // Q30 at 25 kb, SD 3: E[10^(-q/10)] = 1e-3 x exp((3 ln10 / 10)^2 / 2)
    EXPECT_NEAR(errors / 300.0 / 25000, 1.27e-3, 2e-4);
    EXPECT_NEAR(length_sum / 300, 25000, 25);

    // The flow model: a read's bases average its quality (Q_MEAN, SD Q_SD), whatever its length.
    auto const flow = LongReadSetup::Parse("ultima:300:40:25:2");
    double phred_sum = 0, bases = 0;
    errors = 0;
    expected = 0;
    for (int i = 0; i < 3000; ++i) {
        std::string const templ = protal::test::RandomSequence(300, gen);
        hifi::Mutate(templ, rng, flow, read, &errors);
        for (char c : read.qual) {
            phred_sum += Phred(c);
            expected += std::pow(10.0, -Phred(c) / 10.0);
        }
        bases += read.qual.size();
    }
    EXPECT_NEAR(phred_sum / bases, 25.0, 0.6);
    EXPECT_NEAR(errors / expected, 1.0, 0.15);

    // The same stream, the same read.
    LongRng r1(9), r2(9);
    MadeRead m1, m2;
    std::string const templ = protal::test::RandomSequence(5000, gen);
    hifi::Mutate(templ, r1, setup, m1);
    hifi::Mutate(templ, r2, setup, m2);
    EXPECT_EQ(m1.seq, m2.seq);
    EXPECT_EQ(m1.qual, m2.qual);
}

TEST(LongReadSimulation, QshmmModel) {
    // A model of one state per accuracy level 71-99, every base at Q10 (p = 0.1); level 100 has none, and pbsim3 gave
    // its reads Q93 throughout and no error: here it is not drawn, the levels with an HMM keep their weights.
    ScratchDir dir("qshmm");
    std::ostringstream model;
    for (int level = 71; level <= 99; ++level) {
        model << level << " IP 1 1.0\n" << level << " EP 1";
        for (int q = 0; q < 94; ++q) model << ' ' << (q == 10 ? "1.0" : "0");
        model << '\n' << level << " TP 1 1.0\n";
    }
    auto const file = dir.Write("TEST.model", model.str());
    QshmmModel qshmm(file, 0.97, 39, 24, 36);
    EXPECT_EQ(qshmm.AccuracyMin(), 72);   // floor(97 x 0.75)
    EXPECT_EQ(qshmm.AccuracyMax(), 100);  // floor(97 x 1.05), at most 100
    EXPECT_TRUE(qshmm.HasHmm(99));
    EXPECT_FALSE(qshmm.HasHmm(100));

    LongRng rng(4);
    std::map<int, int> levels;
    for (int i = 0; i < 100000; ++i) ++levels[qshmm.DrawAccuracy(rng)];
    EXPECT_EQ(levels.begin()->first, 72);
    EXPECT_EQ(levels.rbegin()->first, 99);  // 100 has no HMM: not drawn
    // weights exp(0.22 x level) over 72-99: the top level takes (1 - e^-0.22) / (1 - e^(-0.22 x 28)) of the reads
    double const top = (1 - std::exp(-0.22)) / (1 - std::exp(-0.22 * 28));
    EXPECT_NEAR(levels[99] / 100000.0, top, 0.01);
    EXPECT_NEAR(static_cast<double>(levels[98]) / levels[99], std::exp(-0.22), 0.03);

    std::mt19937 gen(8);
    MadeRead read;
    std::uint64_t errors = 0, reads = 0, q10_bases = 0, q10_templ = 0;
    for (int i = 0; i < 3000; ++i) {
        std::string const templ = protal::test::RandomSequence(2000, gen);
        qshmm.Mutate(templ, rng, read, &errors);
        ASSERT_EQ(read.seq.size(), read.qual.size());
        ASSERT_EQ(read.qual.find_first_not_of('+'), std::string::npos);  // Q10 throughout: no read at Q93
        ++reads;
        q10_bases += read.seq.size();
        q10_templ += templ.size();
    }
    EXPECT_EQ(reads, 3000u);
    // At Q10: substitutions and insertions 0.1 x (39 + 24) / 99 per emitted base, deletions about 0.1 x 36 / 99 after
    // each: the read's length against its template's about 1 + 0.1 x (24 - 36) / 99.
    EXPECT_NEAR(static_cast<double>(q10_bases) / q10_templ, 1 + 0.1 * (24 - 36) / 99.0, 0.004);
    EXPECT_NEAR(static_cast<double>(errors) / q10_bases, 0.1, 0.01);

    // A model without an HMM in the range: uniform qualities, but never level 100 (Q93 throughout, no error).
    auto const low = dir.Write("LOW.model", "50 IP 1 1.0\n50 EP 1 1.0\n50 TP 1 1.0\n");
    QshmmModel uniform(low, 0.97, 39, 24, 36);
    for (int i = 0; i < 5000; ++i) ASSERT_LT(uniform.DrawAccuracy(rng), 100);
    uniform.Mutate(protal::test::RandomSequence(2000, gen), rng, read);
    EXPECT_NE(read.qual.find_first_not_of('~'), std::string::npos);

    EXPECT_THROW(QshmmModel(dir / "missing.model", 0.97, 39, 24, 36), std::runtime_error);
    // States beyond 50 are written where pbsim3 writes them, but nothing beyond its tables; unknown records fail.
    for (auto const* text : {"100 IP 60 1.0\n", "101 IP 1 1.0\n", "80 XP 1 1.0\n", "80 IP 0 1.0\n"}) {
        auto const broken = dir.Write("BROKEN.model", text);
        EXPECT_THROW(QshmmModel(broken, 0.97, 39, 24, 36), std::runtime_error) << text;
    }
}

// A genome of short contigs gets its weight's share of the bases, as one of a single long contig does: a template is
// placed where it fits (until 2026-10-07 every one was cut at its contig's end, and here the short contigs' genome got
// about a quarter less of each round), and one round is enough.
TEST(LongReadSimulation, ShortContigsGetTheirShare) {
    ScratchDir dir("long share");
    std::mt19937 gen(5);
    std::string fragmented;
    for (int c = 0; c < 30; ++c) fragmented += ">f" + std::to_string(c) + "\n" + protal::test::RandomSequence(2000, gen) + "\n";
    auto const frag = dir.Write("FRAG.fna", fragmented);
    auto const whole = dir.Write("WHOLE.fna", ">w\n" + protal::test::RandomSequence(60000, gen) + "\n");
    LongReadOptions options;
    options.setup = LongReadSetup::Parse("hifi:1000:300:3");
    options.threads = 2;
    std::vector<LongSample> samples(1);
    samples[0] = {"s", dir.path / "s.fq.zst", 2'000'000, 3, {{"FRAG", frag, 60000.0, false}, {"WHOLE", whole, 60000.0, false}}};
    auto const results = SimulateLongReads(samples, options);
    std::uint64_t bases[2] = {0, 0};
    for (auto const& read : ReadFastq(samples[0].out.string())) bases[read.name.rfind("@g1x_", 0) == 0 ? 1 : 0] += read.seq.size();
    EXPECT_NEAR(static_cast<double>(bases[0]) / static_cast<double>(bases[0] + bases[1]), 0.5, 0.03);
    EXPECT_LE(results[0].rounds, 2u);
}

TEST(LongReadSimulation, SamplesTemplatesAndThreads) {
    ScratchDir dir("long reads");
    std::mt19937 gen(21);
    // GA: a contig of 1,500 bases and one of 50 (too short to draw from); GB: 60 kb. A host of two contigs.
    std::string const a1 = protal::test::RandomSequence(1500, gen), a2 = protal::test::RandomSequence(50, gen);
    std::string const b = protal::test::RandomSequence(60000, gen);
    auto const ga = dir.Write("GA.fna", ">a1 x\n" + a1.substr(0, 700) + "\n" + a1.substr(700) + "\n>a2\n" + a2 + "\n");
    auto const gb = dir.Write("GB.fna", ">b\n" + b + "\n");
    std::string const h1 = protal::test::RandomSequence(12000, gen), h2 = protal::test::RandomSequence(8000, gen);
    std::filesystem::create_directories(dir.path / "host");
    dir.Write("host/host.seq", h1 + h2);
    dir.Write("host/host.json", "{\"source\": [\"/x/contigs.fa\", 1, 2], \"contigs\": [[\"chr1\", 0, 12000], "
                                "[\"chr2\", 12000, 8000]], \"bases\": 20000}");

    auto samples_in = [&](std::filesystem::path const& out) {
        std::vector<LongSample> samples(4);
        samples[0] = {"s1", out / "s1.fq.zst", 300000, 7, {{"GA", ga, 1500.0, false}, {"GB", gb, 60000.0, false}}};
        samples[1] = {"s2", out / "s2.fq.gz", 150000, 8, {{"GB", gb, 1.0, false}, {"GA", ga, 1.0, false}}};
        samples[2] = {"s3", out / "s3.fq.zst", 0, 9, {{"GA", ga, 1.0, false}}};
        samples[3] = {"s4", out / "s4.fq.zst", 100000, 10, {{"GA", ga, 1.0, false}, {"host", dir.path / "host", 1.0, true}}};
        return samples;
    };
    LongReadOptions options;
    options.setup = LongReadSetup::Parse("hifi:1000:300:3");  // Q50: about 99% of the reads without an error
    options.threads = 1;
    auto const one = samples_in(dir.path / "one");
    auto const results = SimulateLongReads(one, options);
    options.threads = 4;
    auto const four = samples_in(dir.path / "four");
    auto const results4 = SimulateLongReads(four, options);
    ASSERT_EQ(results.size(), 4u);
    for (std::size_t i = 0; i < one.size(); ++i) {
        // the same files for any number of threads, and no .partial left
        EXPECT_EQ(protal::test::Slurp(one[i].out), protal::test::Slurp(four[i].out)) << one[i].name;
        EXPECT_FALSE(std::filesystem::exists(one[i].out.string() + ".partial"));
        EXPECT_EQ(results[i].reads, results4[i].reads);
    }
    // From a genome store: the same files, when the run writes it and when the next reads it.
    options.genome_store = dir.path / "store";
    for (auto const* run : {"store1", "store2"}) {
        auto const stored = samples_in(dir.path / run);
        SimulateLongReads(stored, options);
        for (std::size_t i = 0; i < one.size(); ++i) {
            EXPECT_EQ(protal::test::Slurp(stored[i].out), protal::test::Slurp(one[i].out)) << run << " " << one[i].name;
        }
    }
    options.genome_store.clear();

    auto contained = [](std::string const& read, std::vector<std::string const*> const& genome) {
        for (auto const* seq : genome) {
            if (seq->find(read) != std::string::npos) return 1;
            if (seq->find(Revcomp(read)) != std::string::npos) return -1;
        }
        return 0;
    };
    for (std::size_t s : {0, 1}) {
        auto const reads = ReadFastq(one[s].out.string());
        ASSERT_EQ(reads.size(), results[s].reads);
        EXPECT_GE(results[s].template_bases, one[s].bases);
        EXPECT_LT(results[s].template_bases, one[s].bases + 20000);
        EXPECT_GE(results[s].rounds, 1u);
        std::map<int, int> per_genome;
        int exact = 0, forward = 0, reverse = 0;
        for (std::size_t r = 0; r < reads.size(); ++r) {
            // g<i>x_<n>: n = 1, 2, ... in the file's order
            auto const& name = reads[r].name;
            auto const x = name.find("x_");
            ASSERT_EQ(name.substr(0, 2), "@g");
            ASSERT_NE(x, std::string::npos);
            EXPECT_EQ(name.substr(x + 2), std::to_string(r + 1));
            int const g = std::stoi(name.substr(2, x - 2));
            ++per_genome[g];
            bool const from_a = one[s].genomes[g].name == "GA";
            if (from_a) EXPECT_LE(reads[r].seq.size(), 1510u) << "a read ends where its contig does";
            int const where = from_a ? contained(reads[r].seq, {&a1}) : contained(reads[r].seq, {&b});
            exact += where != 0;
            forward += where > 0;
            reverse += where < 0;
        }
        EXPECT_GE(exact, 0.9 * reads.size());
        EXPECT_GT(forward, 0.3 * reads.size());
        EXPECT_GT(reverse, 0.3 * reads.size());
        if (s == 0) {  // by weight: GB 60000 of 61500 of the reads' starts
            EXPECT_NEAR(static_cast<double>(per_genome[1]) / reads.size(), 60000.0 / 61500, 0.03);
        }
    }
    EXPECT_EQ(results[2].reads, 0u);
    EXPECT_TRUE(ReadFastq(one[2].out.string()).empty());
    auto const host_reads = ReadFastq(one[3].out.string());
    int from_host = 0, in_host = 0;
    for (auto const& read : host_reads) {
        if (read.name.rfind("@g1x_", 0) != 0) continue;
        ++from_host;
        in_host += contained(read.seq, {&h1, &h2}) != 0;
    }
    EXPECT_GT(from_host, 0.3 * host_reads.size());
    EXPECT_GE(in_host, 0.9 * from_host);

    // Another seed, other reads; a genome without a contig of 100 bases fails the run and leaves no file.
    auto other = samples_in(dir.path / "other");
    other[0].seed = 70;
    options.threads = 2;
    SimulateLongReads({other[0]}, options);
    EXPECT_NE(protal::test::Slurp(other[0].out), protal::test::Slurp(one[0].out));
    auto const tiny = dir.Write("TINY.fna", ">t\n" + a2 + "\n");
    LongSample failing{"bad", dir.path / "bad" / "bad.fq.zst", 10000, 1, {{"TINY", tiny, 1.0, false}}};
    try {
        SimulateLongReads({failing}, options);
        ADD_FAILURE() << "expected a failure";
    } catch (std::runtime_error const& e) {
        EXPECT_NE(std::string(e.what()).find("no sequence of 100 bases"), std::string::npos) << e.what();
    }
    EXPECT_FALSE(std::filesystem::exists(failing.out));
    EXPECT_FALSE(std::filesystem::exists(failing.out.string() + ".partial"));
    options.setup = LongReadSetup::Parse("qshmm:M:8000:6000:0.97");
    EXPECT_THROW(SimulateLongReads({other[0]}, options), std::runtime_error);  // no model file

    // The task files of simulate_metagenomes --long_samples.
    auto const samples_tsv = dir.Write("samples.tsv", "sample\tout\tbases\tseed\nx\t/o/x.fq.zst\t500\t3\ny\t/o/y.fq.gz\t0\t4\n");
    auto const genomes_tsv = dir.Write("genomes.tsv", "sample\tgenome\tfasta\tweight\thost\nx\tGA\t/g/a.fna.gz\t2.5\t0\n"
                                                       "x\thost\t/h\t1\t1\ny\tGB\t/g/b.fna\t1e3\t0\n");
    auto const tasks = ReadLongSamples(samples_tsv, genomes_tsv);
    ASSERT_EQ(tasks.size(), 2u);
    EXPECT_EQ(tasks[0].bases, 500u);
    EXPECT_EQ(tasks[0].seed, 3u);
    ASSERT_EQ(tasks[0].genomes.size(), 2u);
    EXPECT_EQ(tasks[0].genomes[0].weight, 2.5);
    EXPECT_TRUE(tasks[0].genomes[1].host);
    EXPECT_EQ(tasks[1].genomes[0].weight, 1000);
    auto const unknown = dir.Write("unknown.tsv", "sample\tgenome\tfasta\tweight\thost\nz\tGA\t/a\t1\t0\n");
    EXPECT_THROW(ReadLongSamples(samples_tsv, unknown), std::runtime_error);
}
