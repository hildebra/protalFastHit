// The Illumina paired-end reads simulate_metagenomes makes in process (src/RandomForest/IlluminaSimulator.h, on
// ReadPipeline.h): the instrument profiles against their sources, the pairs from their fragments, names, a host's
// pairs, the same files for any number of threads, first reads alone, and samples streamed into named pipes.

#include <gtest/gtest.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <future>
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <poll.h>
#include <sys/stat.h>
#include <unistd.h>

#include "IO/ThreadedGzStream.h"
#include "RandomForest/GenomeStore.h"
#include "RandomForest/IlluminaSimulator.h"
#include "RandomForest/MetagenomeSimulator.h"
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

// Named pipes, each read to its end on a thread of its own; Hold returns once none of them reads any more (a reader
// comes back between reads, every 10 ms at least), Release lets them read on.
class PipeReaders {
public:
    explicit PipeReaders(std::vector<std::filesystem::path> const& paths) : m_bytes(paths.size()) {
        for (std::size_t k = 0; k < paths.size(); ++k) m_threads.emplace_back([this, k, path = paths[k]] { Read(k, path); });
    }
    ~PipeReaders() {
        Release();
        for (auto& thread : m_threads) {
            if (thread.joinable()) thread.join();
        }
    }
    void Hold() {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_hold = true;
        m_cv.wait(lock, [&] { return m_held + m_done == m_threads.size(); });
    }
    void Release() {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_hold = false;
        }
        m_cv.notify_all();
    }
    // Each pipe's bytes, once every writer has closed it.
    std::vector<std::string> Join() {
        for (auto& thread : m_threads) thread.join();
        return m_bytes;
    }

private:
    void Read(std::size_t k, std::filesystem::path const& path) {
        int const fd = ::open(path.c_str(), O_RDONLY);  // once a writer opens it
        std::vector<char> buffer(1 << 16);
        for (bool more = fd >= 0; more;) {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                if (m_hold) {
                    ++m_held;
                    m_cv.notify_all();
                    m_cv.wait(lock, [&] { return !m_hold; });
                    --m_held;
                }
            }
            pollfd ready{fd, POLLIN, 0};
            if (::poll(&ready, 1, 10) <= 0) continue;
            ssize_t const n = ::read(fd, buffer.data(), buffer.size());
            if (n > 0) m_bytes[k].append(buffer.data(), static_cast<std::size_t>(n));
            more = n > 0 || (n < 0 && errno == EINTR);
        }
        if (fd >= 0) ::close(fd);
        std::lock_guard<std::mutex> lock(m_mutex);
        ++m_done;
        m_cv.notify_all();
    }

    std::vector<std::string> m_bytes;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_hold = false;
    std::size_t m_held = 0, m_done = 0;
    std::vector<std::thread> m_threads;
};

// Fills a named pipe that a reader holds open to the last byte, as another writer (without blocking): whole pages while
// it has free ones, then single bytes. -> the bytes written ('#'), none if it was full already.
std::size_t FillPipe(std::filesystem::path const& path) {
    int const fd = ::open(path.c_str(), O_WRONLY | O_NONBLOCK);
    EXPECT_GE(fd, 0) << path << ": " << std::strerror(errno);
    if (fd < 0) return 0;
    std::string const page(4096, '#');
    std::size_t filled = 0;
    for (std::size_t const size : {page.size(), std::size_t{1}}) {
        for (ssize_t n; (n = ::write(fd, page.data(), size)) > 0;) filled += static_cast<std::size_t>(n);
    }
    ::close(fd);
    return filled;
}

// A named pipe opened for writing once a reader has opened it (a minute at most), or -1.
int OpenOnceRead(std::filesystem::path const& path) {
    auto const until = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    for (;;) {
        int const fd = ::open(path.c_str(), O_WRONLY | O_NONBLOCK);
        if (fd >= 0 || errno != ENXIO || std::chrono::steady_clock::now() > until) return fd;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
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
    LongRng rng(3), rng2(4);
    MadeRead r1, r2;
    int close = 0, both_strands = 0;
    for (int i = 0; i < 400; ++i) {
        std::uint32_t const length = model.FragmentLength(rng, genome.size());
        ASSERT_GE(length, 150u);
        std::string const fragment = genome.substr(rng.Below(genome.size() - length + 1), length);
        model.Pair(fragment, 0.0, rng, rng2, r1, &r2);
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

TEST(IlluminaSimulation, GaussianNumbersAndRareEvents) {
    // The ziggurat's numbers: the moments, the shares beyond 1 and 2.5 SD and beyond its tail's start (3.4426).
    LongRng rng(17);
    int const n = 4'000'000;
    double sum = 0, squares = 0, fourth = 0;
    int above1 = 0, above25 = 0, tail = 0;
    for (int i = 0; i < n; ++i) {
        double const z = rng.Gaussian();
        sum += z;
        squares += z * z;
        fourth += z * z * z * z;
        above1 += z > 1;
        above25 += z > 2.5;
        tail += std::fabs(z) > 3.442619855899;
    }
    EXPECT_NEAR(sum / n, 0.0, 0.002);
    EXPECT_NEAR(squares / n, 1.0, 0.003);
    EXPECT_NEAR(fourth / n, 3.0, 0.03);
    EXPECT_NEAR(above1 / static_cast<double>(n), 0.158655, 0.0006);
    EXPECT_NEAR(above25 / static_cast<double>(n), 0.0062097, 0.0002);
    EXPECT_NEAR(tail / static_cast<double>(n), 5.7607e-4, 0.00005);

    // Rare events at exaggerated rates, drawn as geometric gaps: insertions, deletions and no calls per base at their
    // rates; the low state's share at its stationary value (entry 0.05 a cycle, exit 0.25: 1/6, from a high start).
    auto setup = MakeSetup("HS20", 150);
    setup.profile.insertion = setup.profile.deletion = 0.01;
    setup.profile.homopolymer_indels = 1.0;
    setup.profile.n_rate = setup.profile.n_first = 0.01;
    setup.profile.collapse = 0;
    setup.profile.low_rate = setup.profile.low_rate_r2 = 0.05;
    setup.profile.low_growth = 0;
    setup.mean_quality = 40;
    auto const report = ReportProfile(setup, 4000);
    for (int r = 0; r < 2; ++r) {
        EXPECT_NEAR(report.insertions[r], 0.01, 0.001) << r;
        EXPECT_NEAR(report.deletions[r], 0.01, 0.001) << r;
        EXPECT_NEAR(report.ns[r], 0.0099, 0.001) << r;
        double low = 0;  // the low state's bases (N(15, 5): 98% at Q25 or less) and the no calls (Q2)
        for (int q = 0; q <= 25; ++q) low += report.quality_share[r][q];
        double const stationary = 0.05 / 0.30;  // cycle c low with stationary x (1 - 0.7^(c + 1))
        double const mean_low = stationary * (1 - 0.7 * (1 - std::pow(0.70, 150)) / (0.30 * 150));
        EXPECT_NEAR(low, mean_low * 0.982 + report.ns[r], 0.01) << r;
    }
}

// A replay (simulate_metagenomes --from_manifest) of a run's manifest makes the same reads whatever its seed: the
// genomes' art_seed and each sample's run_seed and host_seed are in the manifest; a per-sample manifest replays its
// sample's reads. With first reads only no _R2 is written or named.
TEST(IlluminaSimulation, ReplayFromManifestsMakesTheSameReads) {
    ScratchDir dir("replay");
    std::mt19937 gen(31);
    std::vector<GenomeRecord> genomes;
    for (int g = 0; g < 4; ++g) {
        std::string const fasta = dir.Write("G" + std::to_string(g) + ".fna",
                                            ">c" + std::to_string(g) + "\n" + protal::test::RandomSequence(8000, gen) + "\n");
        genomes.push_back({"G" + std::to_string(g), "d__B;p__P;c__C;o__O;f__F;g__G" + std::to_string(g) + ";s__G" +
                           std::to_string(g) + " sp", fasta, 8000});
    }
    std::filesystem::create_directories(dir.path / "host");
    dir.Write("host/host.seq", protal::test::RandomSequence(20000, gen));
    dir.Write("host/host.json", "{\"source\": [\"/x\", 1, 2], \"contigs\": [[\"chr1\", 0, 20000]], \"bases\": 20000}");
    IlluminaOptions illumina;
    illumina.sequencer = "HSXt";
    illumina.host_folder = dir.path / "host";
    illumina.host_pairs = {40};
    ProfileDesignOptions profile;
    profile.total_read_pairs = 600;
    profile.species_per_sample = 3;
    MetagenomeSimulator original(genomes, illumina, 5);
    auto const samples = original.simulate_samples(profile, 3, "s", dir.path / "original");
    write_combined_manifest(samples, dir.path / "manifest.tsv");
    write_sample_manifest(samples[1], dir.path / "s_2.tsv");

    auto design = read_manifest(dir.path / "manifest.tsv");
    ASSERT_EQ(design.size(), 3u);
    EXPECT_EQ(design[1].run_seed, samples[1].run_seed);
    EXPECT_EQ(design[1].host_seed, samples[1].host_seed);
    MetagenomeSimulator replay(genomes, illumina, 999);  // another seed
    auto const again = replay.replay_samples(design, dir.path / "replay");
    for (std::size_t s = 0; s < 3; ++s) {
        EXPECT_EQ(protal::test::Slurp(again[s].read1_path), protal::test::Slurp(samples[s].read1_path)) << s;
        EXPECT_EQ(protal::test::Slurp(again[s].read2_path), protal::test::Slurp(samples[s].read2_path)) << s;
    }
    MetagenomeSimulator single(genomes, illumina, 1234);  // sample 2 alone: the first of its replay
    auto const one = single.replay_samples(read_manifest(dir.path / "s_2.tsv"), dir.path / "single");
    ASSERT_EQ(one.size(), 1u);
    EXPECT_EQ(protal::test::Slurp(one[0].read1_path), protal::test::Slurp(samples[1].read1_path));
    EXPECT_EQ(protal::test::Slurp(one[0].read2_path), protal::test::Slurp(samples[1].read2_path));

    // first reads only: the same _R1, no _R2 file or name
    illumina.first_reads_only = true;
    MetagenomeSimulator first(genomes, illumina, 5);
    auto const firsts = first.simulate_samples(profile, 3, "s", dir.path / "first");
    for (std::size_t s = 0; s < 3; ++s) {
        EXPECT_TRUE(firsts[s].read2_path.empty());
        EXPECT_EQ(protal::test::Slurp(firsts[s].read1_path), protal::test::Slurp(samples[s].read1_path)) << s;
    }
    EXPECT_FALSE(std::filesystem::exists(dir.path / "first" / "reads" / "s_1_R2.fq.gz"));
    write_combined_manifest(firsts, dir.path / "first.tsv");
    std::istringstream rows(protal::test::Slurp(dir.path / "first.tsv"));
    std::string line;
    std::getline(rows, line);
    while (std::getline(rows, line)) {
        std::vector<std::string> fields;
        std::stringstream split(line);
        for (std::string field; std::getline(split, field, '\t');) fields.push_back(field);
        ASSERT_GE(fields.size(), 10u);
        EXPECT_EQ(fields[9], "") << "fastq_r2 of " << fields[0];
    }
}

// A run that fails while it streams into named pipes ends them, after all it wrote, with what no reader takes for the
// end of a sample: a cut zstd frame (or gzip member), or a FASTQ record without its sequence (plain pipes). Here with
// the pipes full and their readers held when the run fails: the failing genome is a named pipe too, which the run opens
// once GA's pairs are written; the test then holds the readers, fills the pipes up (as another writer) and closes the
// genome's pipe empty. The markers have to wait for the readers to read on.
TEST(IlluminaSimulation, AFailedStreamIsCutOff) {
    ScratchDir dir("failed stream");
    std::mt19937 gen(41);
    auto const ga = dir.Write("GA.fna", ">a\n" + protal::test::RandomSequence(20000, gen) + "\n");
    auto const gone = dir.path / "gone.fna";
    ASSERT_EQ(::mkfifo(gone.c_str(), 0600), 0);
    PairedOptions options;
    options.setup = MakeSetup("NovaSeq", 150);
    options.threads = 1;
    for (bool plain : {false, true}) {
        auto const out = dir.path / (plain ? "plain" : "zstd");
        std::filesystem::create_directories(out);
        std::vector<PairedSample> samples(1);
        // GA's pairs are written, then GONE, without a sequence, fails the run
        samples[0] = {"s", out / "s_R1.fq.zst", out / "s_R2.fq.zst", {{"GA", ga, 4000, 1}, {"GONE", gone, 100, 2}}};
        ASSERT_EQ(::mkfifo(samples[0].r1.c_str(), 0600), 0);
        ASSERT_EQ(::mkfifo(samples[0].r2.c_str(), 0600), 0);
        PipeReaders readers({samples[0].r1, samples[0].r2});
        auto filled = std::async(std::launch::async, [&] {
            int const fd = OpenOnceRead(gone);  // the run reads GONE: GA's pairs are written
            readers.Hold();
            std::vector<std::size_t> padding{FillPipe(samples[0].r1), FillPipe(samples[0].r2)};
            if (fd >= 0) ::close(fd);  // the run fails, its pipes full
            std::this_thread::sleep_for(std::chrono::milliseconds(200));  // time for a writer to give up on them
            readers.Release();
            return padding;
        });
        options.plain_pipes = plain;
        EXPECT_THROW(SimulatePairs(samples, options), std::runtime_error);
        ASSERT_EQ(filled.wait_for(std::chrono::seconds(60)), std::future_status::ready);
        auto const padding = filled.get();
        auto const got = readers.Join();
        std::string const marker = plain ? std::string("@simulate_metagenomes_failed\n")
                                         : std::string{'\x28', '\xb5', '\x2f', '\xfd', '\x24'};
        for (std::size_t r = 0; r < 2; ++r) {
            std::string const tail = std::string(padding[r], '#') + marker;
            ASSERT_TRUE(got[r].ends_with(tail)) << r << ": the marker, after all that was in the pipe";
            std::string const pairs = got[r].substr(0, got[r].size() - tail.size());  // written before the failure
            std::string const text = plain ? pairs : ReadWholeFile(dir.Write("pairs.fq.zst", pairs));
            EXPECT_EQ(std::count(text.begin(), text.end(), '\n'), 4 * 4000) << r << ": all of GA's pairs";
            if (!plain) {
                EXPECT_THROW(ReadWholeFile(dir.Write("cut.fq.zst", pairs + marker)), std::runtime_error) << "a cut zstd frame";
            }
        }
    }
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
    // From a genome store: the same files, when the run writes it and when the next reads it.
    options.genome_store = dir.path / "store";
    for (auto const* run : {"store1", "store2"}) {
        auto const stored = samples_in(dir.path / run, ".fq.zst");
        SimulatePairs(stored, options);
        for (std::size_t s = 0; s < 3; ++s) {
            EXPECT_EQ(protal::test::Slurp(stored[s].r1), protal::test::Slurp(one[s].r1)) << run;
            if (!one[s].r2.empty()) EXPECT_EQ(protal::test::Slurp(stored[s].r2), protal::test::Slurp(one[s].r2)) << run;
        }
    }
    options.genome_store.clear();

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

    // With plain_pipes: plain FASTQ through the pipes, the text the files hold; a regular file stays compressed.
    std::filesystem::create_directories(dir.path / "plain");
    auto plain = samples_in(dir.path / "plain", ".fq.zst");
    plain.resize(1);
    ASSERT_EQ(::mkfifo(plain[0].r1.c_str(), 0600), 0);
    ASSERT_EQ(::mkfifo(plain[0].r2.c_str(), 0600), 0);
    auto drain = [](std::filesystem::path path) {  // each pipe on a thread of its own: no lockstep needed
        return std::async(std::launch::async, [path] { return protal::test::Slurp(path); });
    };
    auto text1 = drain(plain[0].r1), text2 = drain(plain[0].r2);
    options.plain_pipes = true;
    SimulatePairs(plain, options);
    ASSERT_EQ(text1.wait_for(std::chrono::seconds(60)), std::future_status::ready);
    ASSERT_EQ(text2.wait_for(std::chrono::seconds(60)), std::future_status::ready);
    EXPECT_EQ(text1.get(), ReadWholeFile(one[0].r1));
    EXPECT_EQ(text2.get(), ReadWholeFile(one[0].r2));
    auto regular = samples_in(dir.path / "regular", ".fq.zst");
    regular.resize(1);
    SimulatePairs(regular, options);
    EXPECT_EQ(protal::test::Slurp(regular[0].r1), protal::test::Slurp(one[0].r1));
}
