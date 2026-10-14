#!/bin/bash
# Moves the bodies of Options.h's parsing and checks into src/Options.cpp: CxxOptions, ToString, PrintHelp, LoadFromMap,
# PrepareAndCheckValidity and OptionsFromArguments. The class keeps their declarations (and default arguments); every
# body moves verbatim, four spaces less indented. The functions are found by their signatures (usage:
# make_options2.sh COMMIT), and each must end where expected.
set -euo pipefail
COMMIT=${1:?usage: make_options2.sh COMMIT}
W=C:/Users/hildebra/Documents/locDev/protal-split
S=$(cd "$(dirname "$0")" && pwd)
O=$S/Options.upstream.h
git -C $W show $COMMIT:src/Options.h > $O
n=$(wc -l < $O)
r() { sed -n "$1,$2p" $O; }
dedent() { sed 's/^    //'; }
first() {  # first line matching the (extended) regex at or after line $2
    awk -v re="$1" -v from="${2:-1}" 'NR >= from && $0 ~ re { print NR; exit }' $O
}
one() {  # the only line matching the regex
    local c
    c=$(grep -cE -- "$1" $O || true)
    [[ $c == 1 ]] || { echo "'$1' found $c times" >&2; exit 1; }
    grep -nE -- "$1" $O | cut -d: -f1
}
inc=$(one '^#include <cxxopts.hpp>$')
cx=$(one '^    static cxxopts::Options CxxOptions\(\) \{$');                        cx_end=$(first '^    }$' $cx)
ts=$(one '^        std::string ToString\(\) const \{$');                            ts_end=$(first '^        }$' $ts)
ph=$(one '^        void PrintHelp\(bool show_dev=false\) \{$');                      ph_end=$(first '^        }$' $ph)
lm=$(one '^        static bool LoadFromMap\(std::string map_path,');                 lm_sig=$(first '= nullptr\) \{$' $lm); lm_end=$(first '^        }$' $lm_sig)
pc=$(one '^        bool PrepareAndCheckValidity\(bool force_read_check=false\) \{$'); pc_end=$(first '^        }$' $pc)
of=$(one '^        static Options OptionsFromArguments\(int argc, char \*argv\[\]\) \{$'); of_end=$(first '^        }$' $of)
echo "cxxopts include $inc; CxxOptions $cx-$cx_end; ToString $ts-$ts_end; PrintHelp $ph-$ph_end; LoadFromMap $lm-$lm_sig-$lm_end; PrepareAndCheckValidity $pc-$pc_end; OptionsFromArguments $of-$of_end of $n"
# Each in order, and the parser ends the class.
[[ $inc -lt $cx && $cx_end -lt $ts && $ts_end -lt $ph && $ph_end -lt $lm && $lm_end -lt $pc && $pc_end -lt $of ]] || { echo "unexpected order" >&2; exit 1; }
[[ $(sed -n "$((of_end + 1))p" $O) == '    };' ]] || { echo "OptionsFromArguments does not end the class" >&2; exit 1; }
[[ $((lm_sig - lm)) == 9 ]] || { echo "LoadFromMap's signature is not 10 lines" >&2; exit 1; }
# cxxopts only in the moved functions
moved_cx=$( { r $inc $inc; r $cx $cx_end; r $ph $ph_end; r $of $of_end; } | grep -c cxxopts)
[[ $(grep -c cxxopts $O) == $moved_cx ]] || { echo "cxxopts used in Options.h outside the moved functions" >&2; exit 1; }

{
    r 1 $((inc - 1))
    r $((inc + 1)) $((cx - 1))
    echo '    // The options'"'"' definitions (CxxOptions), the parsing and the checks before a run: Options.cpp.'
    r $((cx_end + 1)) $((ts - 1))
    echo '        std::string ToString() const;'
    r $((ts_end + 1)) $((ph - 1))
    echo '        void PrintHelp(bool show_dev=false);'
    r $((ph_end + 1)) $((lm - 1))
    r $lm $lm_sig | sed 's/ = nullptr) {$/ = nullptr);/'
    r $((lm_end + 1)) $((pc - 1))
    echo '        bool PrepareAndCheckValidity(bool force_read_check=false);'
    r $((pc_end + 1)) $((of - 1))
    echo '        static Options OptionsFromArguments(int argc, char *argv[]);'
    r $((of_end + 1)) $n
} > $W/src/Options.h

{
    echo '// The options'"'"' definitions on the command line (CxxOptions), the map file, the parsing and the checks before a run'
    echo '// (Options.h declares them). Here rather than in Options.h, which nearly every unit of protal includes'
    echo '// (docs/claude/2026-10-14-compile-time).'
    echo
    echo '#include <cxxopts.hpp>'
    echo '#include "Options.h"'
    echo
    echo 'namespace protal {'
    r $cx $cx_end
    echo
    echo '    std::string Options::ToString() const {'
    r $((ts + 1)) $ts_end | dedent
    echo
    echo '    void Options::PrintHelp(bool show_dev) {'
    r $((ph + 1)) $ph_end | dedent
    echo
    echo '    bool Options::LoadFromMap(std::string map_path, std::string& output_dir, std::string& strain_output_dir, std::string& misc_output_dir,'
    echo '                              std::vector<std::string>& prefix_list,'
    echo '                              std::vector<std::string>& first_list, std::vector<std::string>& second_list,'
    echo '                              std::vector<std::string>& sam_list, std::vector<std::string>& profile_list,'
    echo '                              std::vector<std::string>& samplenames_list, std::vector<std::string>& profile_truth_list,'
    echo '                              std::vector<std::string>& read_type_list, std::vector<std::string>& unmapped_reads_list,'
    echo '                              std::string const& sam_ending, bool* sams_named,'
    echo '                              std::vector<std::string>* sparse_first_list,'
    echo '                              std::vector<std::string>* flex_indels_list,'
    echo '                              std::vector<std::string>* sparse_first_pick_list) {'
    r $((lm_sig + 1)) $lm_end | dedent
    echo
    echo '    bool Options::PrepareAndCheckValidity(bool force_read_check) {'
    r $((pc + 1)) $pc_end | dedent
    echo
    echo '    Options Options::OptionsFromArguments(int argc, char *argv[]) {'
    r $((of + 1)) $of_end | dedent
    echo '}'
} > $W/src/Options.cpp
# The definition's parameters are the declaration's, less their defaults.
diff <(r $lm $lm_sig | sed -E 's/^ +//; s/static //; s/ = (kDefaultSamEnding|nullptr)//g; s/ \{$//' | tr -s ' ' | tr -d '\n') \
     <(grep -A9 '^    bool Options::LoadFromMap(' $W/src/Options.cpp | sed -E 's/^ +//; s/Options:://; s/ \{$//' | tr -s ' ' | tr -d '\n') \
    && echo "LoadFromMap's parameters match"
echo "Options.h $n -> $(wc -l < $W/src/Options.h); Options.cpp $(wc -l < $W/src/Options.cpp)"
