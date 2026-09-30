# Adds a gene-page log to a COPY of the source tree (never the checkout): with PROTAL_TOUCH=1 the
# build counts the Sequence() calls and the 4 KB pages of the gene arena they hand out, and prints
# both after aligning and at the end. Measurement only.
#   rsync -a --exclude /.git --exclude '/build*' CHECKOUT/ COPY/ && perl touch_patch.pl COPY
#   cmake -S COPY -B build-touch -G Ninja -DCMAKE_BUILD_TYPE=Release && ninja -C build-touch protal
#   PROTAL_TOUCH=1 build-touch/protal ... 2>&1 | grep TOUCH        (use a fresh -o: an existing SAM skips aligning)
use strict; use warnings;
my $root = shift or die "usage: perl touch_patch.pl SOURCE_COPY\n";
sub slurp { open my $in, "<", $_[0] or die "$_[0]: $!"; local $/; my $s = <$in>; close $in; $s }
sub spit { open my $out, ">", $_[0] or die; print $out $_[1]; close $out }

my $f = "$root/src/SequenceUtils/GenomeLoader.h";
my $s = slurp($f);
my $log = <<'LOG';
    // MEASUREMENT ONLY: which 4 KB pages of the gene arena Gene::Sequence() handed out (PROTAL_TOUCH=1).
    struct TouchLog {
        static inline const char* base = nullptr;
        static inline uint64_t bytes = 0;
        static inline std::atomic<uint8_t>* pages = nullptr;
        static inline bool on = false;
        static inline std::atomic<uint64_t> calls{0};
        static void Init(const char* b, uint64_t n) {
            if (!getenv("PROTAL_TOUCH")) return;
            base = b; bytes = n; pages = new std::atomic<uint8_t>[(n >> 12) + 2]();
            on = true;
        }
        static void Mark(const char* p, size_t len) {
            calls.fetch_add(1, std::memory_order_relaxed);
            if (!on || p < base || p >= base + bytes || len == 0) return;
            uint64_t a = (p - base) >> 12, z = (p - base + len - 1) >> 12;
            for (; a <= z; a++) pages[a].store(1, std::memory_order_relaxed);
        }
        static void Report() {
            if (!on) return;
            uint64_t n = 0, t = (bytes >> 12) + 1;
            for (uint64_t i = 0; i < t; i++) n += pages[i].load();
            fprintf(stderr, "Sequence() calls %lu; ", calls.load());
            fprintf(stderr, "TOUCH pages %lu of %lu (%.1f%%), %.1f of %.1f MB\n", n, t, 100.0 * n / t, n * 4096 / 1e6, t * 4096 / 1e6);
        }
    };

LOG
$s =~ s/(    class Gene \{\n        static const size_t DEFAULT)/$log$1/ or die "no class Gene";
$s =~ s/(std::string_view Sequence\(\) const \{\n)(\s+return m_arena \?)/$1            if (m_arena) TouchLog::Mark(m_arena, m_length);\n$2/ or die "no Sequence";
$s =~ s/(                m_arenas\.emplace_back\(std::move\(arena\)\);)/$1\n                TouchLog::Init(m_arenas.back().get(), total);/ or die "no arenas";
$s =~ s/#include <sparse_set.h>/#include <sparse_set.h>\n#include <atomic>/;
spit($f, $s);

$f = "$root/src/RunProtal.h";
$s = slurp($f);
$s =~ s/(            bm_classify\.Stop\(\);\n            bm_classify\.PrintResults\(\);)/$1\n            std::cerr << "after aligning: "; TouchLog::Report();/ or die "no classify";
$s =~ s/(        Benchmark bm_total\("Run protal"\);)/$1\n        struct AtEnd { ~AtEnd() { std::cerr << "at end: "; TouchLog::Report(); } } at_end;/ or die "no total";
spit($f, $s);
print "patched $root\n";
