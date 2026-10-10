// ColumnWeightsMsa.h - column_weights.tsv from GTDB's own alignments and tree (since 2026-10-10): no alignment of
// protal's, each family's columns those of the marker's alignment, and each genus's and family's ancestral sequence
// by parsimony on the species tree.
//
// GTDB distributes, per marker, the representatives' proteins aligned by GTDB-Tk's profile HMM (insertions trimmed:
// every row is its protein less the residues of the model's insert states, on the model's columns), and the tree of
// all representatives it infers from them. gtdb_to_protal_db.py writes them into the database folder: column_msa/
// <geneid>.faa[.zst] (">taxid" and the aligned row) and species_tree.nwk (leaves named by taxid). Since the columns
// are the same for every copy of a gene, nothing has to be aligned: a copy's residues are placed on its row
// (PlaceResidues), each placed residue's codon on the family's three base columns of that amino-acid column. The
// protein alignments of ColumnWeightsBuild.h agreed with these on 98.8% of the residue pairs of a real family, but
// lost the genera beyond their limits; these lose none (docs/claude/2026-10-10-real-ancestry).
//
// Per family and gene (BuildFamilyFromMsa):
//  - within: each genus's species' bases per column (every species of the genus, no reference needed), the genus's
//    raw share where two or more were compared, averaged over the family's genera with the pseudocounts over all the
//    species compared (CodeOfPooled), as ColumnWeightsBuild.h does;
//  - the genera's ancestral sequences: Fitch parsimony on the tree pruned to the genus's species (a genus of one
//    species: its own copy), per base column and per amino-acid column; a tie at the root broken by the states'
//    counts among the species, then by the lowest code;
//  - among: the genus ancestors' bases per column, CodeOf(agree, compared) over every genus of the family (no cap of
//    Settings::genera: the ancestors cost no alignment, and ten genus references capped the code at 5);
//  - aa: the within share of the amino acids and the genus ancestors', pooled;
//  - the consensus base, which polarises a small genus's ancestry sites: the family's ancestral base (Fitch over all
//    its species), where the root's state is unambiguous and kMinConsensusGenera genera or more carry a base.
#pragma once

#include "ColumnWeightsBuild.h"
#include "Zstd.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <istream>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace protal::column_weights::msa {
    inline const std::string kFolder = "column_msa";          // per gene <geneid>.faa or .faa.zst
    inline const std::string kTreeFile = "species_tree.nwk";  // the species tree, leaves named by taxid
    inline constexpr size_t kAnchor = 4;           // residues that must match in a row to place one (PlaceResidues)
    inline constexpr size_t kSearch = 64;          // residues searched ahead in the protein for the next placement
    inline constexpr size_t kFarAnchor = 8;        // residues that must match beyond kSearch (the rest of the protein)
    inline constexpr double kMinPlaced = 0.5;      // of a row's residues placed, for the copy to be used
    inline constexpr uint32_t kBaseStates = 4, kAminoStates = 21;  // A, C, G, T; the 20 amino acids and the stop

    // The species tree: per node its parent (-1 at the root) and children; the leaves by taxid.
    struct Tree {
        std::vector<int32_t> parent;
        std::vector<std::vector<int32_t>> children;
        std::unordered_map<uint32_t, int32_t> leaf_of;

        bool Empty() const { return parent.empty(); }

        // Newick: labels of leaves that are numbers become taxids; other labels, quoted labels, branch lengths and
        // support values are skipped. Empty on a malformed tree, with `error` set.
        static Tree Parse(std::string const& text, std::string& error) {
            Tree t;
            std::vector<int32_t> stack;
            size_t i = 0;
            auto new_node = [&](int32_t up) {
                t.parent.push_back(up);
                t.children.emplace_back();
                if (up >= 0) t.children[static_cast<size_t>(up)].push_back(static_cast<int32_t>(t.parent.size() - 1));
                return static_cast<int32_t>(t.parent.size() - 1);
            };
            auto skip_label = [&](int32_t node) {
                std::string label;
                if (i < text.size() && text[i] == '\'') {
                    size_t const end = text.find('\'', i + 1);
                    i = end == std::string::npos ? text.size() : end + 1;
                } else {
                    while (i < text.size() && text[i] != ',' && text[i] != '(' && text[i] != ')' && text[i] != ':' &&
                           text[i] != ';' && text[i] != ' ' && text[i] != '\n' && text[i] != '\r') {
                        label.push_back(text[i++]);
                    }
                }
                if (i < text.size() && text[i] == ':') {
                    i++;
                    while (i < text.size() && text[i] != ',' && text[i] != ')' && text[i] != ';') i++;
                }
                if (node >= 0 && t.children[static_cast<size_t>(node)].empty() && !label.empty() &&
                    std::all_of(label.begin(), label.end(), [](char c) { return c >= '0' && c <= '9'; })) {
                    t.leaf_of[static_cast<uint32_t>(std::stoul(label))] = node;
                }
            };
            int32_t current = -1;
            while (i < text.size()) {
                char const c = text[i];
                if (c == '(') {
                    current = new_node(stack.empty() ? -1 : stack.back());
                    stack.push_back(current);
                    i++;
                } else if (c == ',') {
                    i++;
                } else if (c == ')') {
                    if (stack.empty()) {
                        error = "unbalanced parentheses";
                        return {};
                    }
                    current = stack.back();
                    stack.pop_back();
                    i++;
                    skip_label(current);
                } else if (c == ';') {
                    break;
                } else if (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
                    i++;
                } else {
                    int32_t const leaf = new_node(stack.empty() ? -1 : stack.back());
                    skip_label(leaf);
                }
            }
            if (!stack.empty() || t.parent.empty()) {
                error = t.parent.empty() ? "no nodes" : "unbalanced parentheses";
                return {};
            }
            return t;
        }

        static Tree Read(std::string const& path, std::string& error) {
            zstd::InputFile in(zstd::Resolve(path));
            if (!in.IsOpen()) {
                error = "cannot open " + path;
                return {};
            }
            std::string text((std::istreambuf_iterator<char>(in.Stream())), std::istreambuf_iterator<char>());
            return Parse(text, error);
        }
    };

    // A tree induced on some leaves (unary nodes contracted): the nodes in post-order (the root last), each with its
    // children (indices into this tree) or, for a leaf, the leaf's index in the taxids it was induced on.
    struct Induced {
        std::vector<std::vector<int32_t>> children;
        std::vector<int32_t> leaf;  // per node: the index into the taxids, -1 for an inner node

        size_t Leaves() const { return static_cast<size_t>(std::count_if(leaf.begin(), leaf.end(), [](int32_t l) { return l >= 0; })); }
    };

    // The tree induced on the taxids that are leaves of `tree` (the others are left out). Costs the leaves' paths to
    // the root, not the tree.
    inline Induced Induce(Tree const& tree, std::vector<uint32_t> const& taxids) {
        Induced out;
        std::unordered_map<int32_t, int32_t> index_of_leaf;   // tree node -> index into taxids
        std::unordered_map<int32_t, std::vector<int32_t>> kids;  // tree node -> its children on the leaves' paths
        std::unordered_set<int32_t> seen;
        int32_t top = -1;
        for (size_t k = 0; k < taxids.size(); k++) {
            auto const it = tree.leaf_of.find(taxids[k]);
            if (it == tree.leaf_of.end()) continue;
            index_of_leaf[it->second] = static_cast<int32_t>(k);
            int32_t node = it->second;
            seen.insert(node);
            while (tree.parent[static_cast<size_t>(node)] >= 0) {
                int32_t const up = tree.parent[static_cast<size_t>(node)];
                bool const fresh = seen.insert(up).second;
                kids[up].push_back(node);
                node = up;
                if (!fresh) break;
            }
            if (tree.parent[static_cast<size_t>(node)] < 0) top = node;
        }
        if (index_of_leaf.empty()) return out;
        // Down from the tree's root (or the only leaf) while a node has one child on the paths.
        int32_t root = top < 0 ? index_of_leaf.begin()->first : top;
        while (!index_of_leaf.count(root) && kids.count(root) && kids[root].size() == 1) root = kids[root][0];
        // Post-order, contracting the nodes with one child on the paths.
        auto contract = [&](int32_t node) {
            while (!index_of_leaf.count(node) && kids.count(node) && kids[node].size() == 1) node = kids[node][0];
            return node;
        };
        std::vector<std::pair<int32_t, size_t>> stack{ { contract(root), 0 } };
        std::vector<std::vector<int32_t>> pending{ {} };  // per stack level the finished children's indices
        while (!stack.empty()) {
            auto& [node, next] = stack.back();
            auto const& ks = kids.count(node) && !index_of_leaf.count(node) ? kids[node] : std::vector<int32_t>{};
            if (next < ks.size()) {
                int32_t const child = contract(ks[next++]);
                stack.emplace_back(child, 0);
                pending.emplace_back();
                continue;
            }
            out.children.push_back(std::move(pending.back()));
            auto const leaf = index_of_leaf.find(node);
            out.leaf.push_back(leaf == index_of_leaf.end() ? -1 : leaf->second);
            pending.pop_back();
            stack.pop_back();
            if (!pending.empty()) pending.back().push_back(static_cast<int32_t>(out.leaf.size() - 1));
        }
        return out;
    }

    // Fitch parsimony per column: `states[k][c]` is leaf k's state at column c (>= `n_states`: unknown, which
    // constrains nothing). The root's state per column (`n_states` where no leaf knows it), a tie broken by the
    // states' counts among the leaves, then by the lowest state; `ambiguous` counts the columns whose root set held
    // more than one state.
    inline std::vector<uint8_t> Fitch(Induced const& tree, std::vector<std::vector<uint8_t>> const& states, size_t columns,
                                      uint32_t n_states, size_t& ambiguous) {
        uint32_t const all = (n_states >= 32 ? ~0u : (1u << n_states) - 1);
        std::vector<std::vector<uint32_t>> mask(tree.leaf.size());
        for (size_t n = 0; n < tree.leaf.size(); n++) {
            auto& m = mask[n];
            m.assign(columns, all);
            if (tree.leaf[n] >= 0) {
                auto const& s = states[static_cast<size_t>(tree.leaf[n])];
                for (size_t c = 0; c < columns && c < s.size(); c++) {
                    if (s[c] < n_states) m[c] = 1u << s[c];
                }
                continue;
            }
            bool first = true;
            for (int32_t const child : tree.children[n]) {
                auto const& cm = mask[static_cast<size_t>(child)];
                for (size_t c = 0; c < columns; c++) {
                    if (first) {
                        m[c] = cm[c];
                    } else {
                        uint32_t const both = m[c] & cm[c];
                        m[c] = both ? both : (m[c] | cm[c]);
                    }
                }
                first = false;
            }
            for (int32_t const child : tree.children[n]) std::vector<uint32_t>().swap(mask[static_cast<size_t>(child)]);
        }
        std::vector<uint8_t> root(columns, static_cast<uint8_t>(n_states));
        if (mask.empty()) return root;
        auto const& rm = mask.back();
        for (size_t c = 0; c < columns; c++) {
            if (rm[c] == all) {
                bool known = false;
                for (auto const& s : states) known = known || (c < s.size() && s[c] < n_states);
                if (!known) continue;
            }
            std::array<uint32_t, 32> counts{};
            for (auto const& s : states) {
                if (c < s.size() && s[c] < n_states) counts[s[c]]++;
            }
            uint8_t best = static_cast<uint8_t>(n_states);
            size_t in_set = 0;
            for (uint32_t s = 0; s < n_states; s++) {
                if (!(rm[c] >> s & 1u)) continue;
                in_set++;
                if (best == n_states || counts[s] > counts[best]) best = static_cast<uint8_t>(s);
            }
            ambiguous += in_set > 1;
            root[c] = best;
        }
        return root;
    }

    // Per residue of `protein` the column of GTDB's aligned `row` it sits on (-1: not on a column, an insertion the
    // row trimmed or a residue the row lacks). The row is the protein less its insertions, so its residues are found
    // in order: each is placed where it and the next kAnchor - 1 residues of the row match the protein (fewer at the
    // row's end), searched up to kSearch residues ahead, then through the rest of the protein with kFarAnchor
    // residues (past a long insertion or N-terminal extension); where no window matches (a residue just before an
    // insertion, whose next residues follow the insertion in the protein), right after the previous placed residue if
    // it matches there; a row residue found nowhere (a start codon GTDB translated otherwise) is skipped. Returns how
    // many of the row's residues were placed.
    inline size_t PlaceResidues(std::string const& protein, std::string const& row, std::vector<int32_t>& column_of) {
        column_of.assign(protein.size(), -1);
        std::vector<std::pair<char, int32_t>> residues;
        for (size_t c = 0; c < row.size(); c++) {
            char const a = row[c];
            if (a == '-' || a == '.') continue;
            residues.emplace_back(a >= 'a' && a <= 'z' ? static_cast<char>(a - 32) : a, static_cast<int32_t>(c));
        }
        size_t placed = 0, j = 0;
        bool previous = false;  // the previous row residue was placed (at j - 1)
        for (size_t i = 0; i < residues.size() && j < protein.size(); i++) {
            // Near (kSearch residues ahead) with kAnchor residues; beyond, past a long insertion or N-terminal extension,
            // with kFarAnchor, against chance matches.
            auto search = [&](size_t from, size_t to, size_t anchor) {
                size_t const want = std::min(anchor, residues.size() - i);
                for (size_t k = from; k < protein.size() && k < to; k++) {
                    if (k + want > protein.size()) break;
                    bool match = true;
                    for (size_t t = 0; t < want && match; t++) match = protein[k + t] == residues[i + t].first;
                    if (match) return k;
                }
                return protein.size();
            };
            size_t found = search(j, j + kSearch, kAnchor);
            // Beyond, the longer window only (a 4-residue match there, even a unique one, was a chance hit in a
            // low-complexity stretch often enough to cost 1.4% of a real family's residue pairs).
            if (found == protein.size() && j + kSearch < protein.size()) found = search(j + kSearch, protein.size(), kFarAnchor);
            if (found == protein.size() && previous && protein[j] == residues[i].first) found = j;
            previous = found != protein.size();
            if (!previous) continue;  // this row residue is not in the protein here: skipped
            column_of[found] = residues[i].second;
            placed++;
            j = found + 1;
        }
        return placed;
    }

    // The aligned rows of one gene: column_msa/<gene>.faa[.zst] in `folder`, the rows of `wanted` taxids. Empty when
    // the file is missing; `width` gets the rows' length.
    inline std::unordered_map<uint32_t, std::string> ReadRows(std::string const& folder, uint32_t gene,
                                                              std::unordered_set<uint32_t> const& wanted, size_t& width) {
        std::unordered_map<uint32_t, std::string> rows;
        width = 0;
        std::string const path = zstd::Resolve(folder + "/" + std::to_string(gene) + ".faa");
        if (!std::filesystem::exists(path)) return rows;
        zstd::InputFile in(path);
        std::string line;
        uint32_t taxid = 0;
        bool keep = false;
        while (std::getline(in.Stream(), line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (line[0] == '>') {
                taxid = static_cast<uint32_t>(std::strtoul(line.c_str() + 1, nullptr, 10));
                keep = wanted.count(taxid) > 0;
                if (keep) rows[taxid].clear();
            } else if (keep) {
                rows[taxid] += line;
            }
        }
        for (auto const& [t, r] : rows) width = std::max(width, r.size());
        return rows;
    }

    struct MsaStats {
        size_t copies = 0;          // copies of a family's species considered
        size_t without_row = 0;     // with no row in GTDB's alignment
        size_t badly_placed = 0;    // whose row's residues were placed below kMinPlaced
        size_t residues = 0, placed = 0;  // the used copies' residues, those on a column
        size_t genus_ancestors = 0, genus_trees = 0;  // genera with an ancestor; of them from a tree (2+ species on it)
        size_t family_columns = 0, family_ambiguous = 0;  // the families' ancestral base columns, ambiguous at the root
        size_t polarised = 0;       // columns with a consensus (the family's ancestral) base

        void Add(MsaStats const& o) {
            copies += o.copies; without_row += o.without_row; badly_placed += o.badly_placed;
            residues += o.residues; placed += o.placed;
            genus_ancestors += o.genus_ancestors; genus_trees += o.genus_trees;
            family_columns += o.family_columns; family_ambiguous += o.family_ambiguous; polarised += o.polarised;
        }
    };

    // One family's copies of a gene (as BuildFamily takes them) and their rows (`rows[i]`: copy i's aligned row,
    // empty if none; all of one width): the family row from GTDB's columns and every placed copy's mapping onto them.
    inline bool BuildFamilyFromMsa(uint32_t family, uint32_t gene, std::vector<uint32_t> const& taxids,
                                   std::vector<std::string> const& seqs, std::vector<uint32_t> const& genus_of,
                                   std::vector<std::string const*> const& rows, size_t width, Tree const& tree,
                                   MsaStats& stats, std::vector<Family>& families,
                                   std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>>& copies) {
        if (width == 0 || 3 * width > kMaxLength) return false;
        size_t const length = 3 * width;
        // Each copy's bases on the base columns and amino acids on the amino-acid columns (4, kAminoAcids - 1: none).
        std::vector<size_t> used;
        std::vector<std::vector<uint8_t>> bases(taxids.size()), aminos(taxids.size());
        std::vector<std::vector<Run>> runs(taxids.size());
        std::vector<int32_t> column_of;
        for (size_t i = 0; i < taxids.size(); i++) {
            if (genus_of[i] == 0 || seqs[i].empty()) continue;
            stats.copies++;
            if (!rows[i] || rows[i]->empty()) {
                stats.without_row++;
                continue;
            }
            std::string const protein = Translate(seqs[i]);
            size_t const row_residues = static_cast<size_t>(std::count_if(rows[i]->begin(), rows[i]->end(),
                                                                         [](char c) { return c != '-' && c != '.'; }));
            size_t const placed = PlaceResidues(protein, *rows[i], column_of);
            if (row_residues == 0 || static_cast<double>(placed) < kMinPlaced * static_cast<double>(row_residues)) {
                stats.badly_placed++;
                continue;
            }
            stats.residues += protein.size();
            stats.placed += placed;
            bases[i].assign(length, kNoBase);
            aminos[i].assign(width, static_cast<uint8_t>(kAminoAcids - 1));
            std::vector<int32_t> copy_to_column(seqs[i].size(), -1);
            for (size_t r = 0; r < column_of.size(); r++) {
                if (column_of[r] < 0) continue;
                size_t const col = static_cast<size_t>(column_of[r]);
                if (col >= width) continue;
                for (size_t k = 0; k < 3 && 3 * r + k < seqs[i].size(); k++) {
                    bases[i][3 * col + k] = BaseCode(seqs[i][3 * r + k]);
                    copy_to_column[3 * r + k] = static_cast<int32_t>(3 * col + k);
                }
                aminos[i][col] = detail::AminoIndex(protein[r]);
            }
            runs[i] = RunsOf(copy_to_column);
            used.push_back(i);
        }
        if (used.empty()) return false;

        std::map<uint32_t, std::vector<size_t>> genera;
        for (size_t const i : used) genera[genus_of[i]].push_back(i);
        // within: the genera's raw shares, pooled; the genus ancestors for among and aa.
        std::vector<double> within_sum(length, 0.0), aa_sum(width, 0.0);
        std::vector<uint32_t> within_n(length, 0), within_compared(length, 0), aa_n(width, 0), aa_compared(width, 0);
        std::vector<std::array<uint32_t, 4>> among_votes(length, std::array<uint32_t, 4>{});
        std::vector<std::array<uint32_t, kAminoAcids>> among_aa(width, std::array<uint32_t, kAminoAcids>{});
        std::vector<uint32_t> genera_with_base(length, 0);
        for (auto const& [genus, members] : genera) {
            for (size_t c = 0; c < length; c++) {
                std::array<uint32_t, 4> v{};
                for (size_t const i : members) if (bases[i][c] < 4) v[bases[i][c]]++;
                uint32_t const n = v[0] + v[1] + v[2] + v[3];
                genera_with_base[c] += n > 0;
                if (n < kMinCompared) continue;
                within_sum[c] += static_cast<double>(detail::Best(v).first) / static_cast<double>(n);
                within_n[c]++;
                within_compared[c] += n;
            }
            for (size_t c = 0; c < width; c++) {
                std::array<uint32_t, kAminoAcids - 1> v{};
                uint32_t n = 0;
                for (size_t const i : members) {
                    if (aminos[i][c] < kAminoAcids - 1) {
                        v[aminos[i][c]]++;
                        n++;
                    }
                }
                if (n < kMinCompared) continue;
                aa_sum[c] += static_cast<double>(detail::Best(v).first) / static_cast<double>(n);
                aa_n[c]++;
                aa_compared[c] += n;
            }
            // The genus's ancestral sequence: its own copy for one species, else Fitch on the tree pruned to it.
            std::vector<uint32_t> member_taxids;
            std::vector<std::vector<uint8_t>> member_bases, member_aminos;
            for (size_t const i : members) {
                member_taxids.push_back(taxids[i]);
                member_bases.push_back(bases[i]);
                member_aminos.push_back(aminos[i]);
            }
            std::vector<uint8_t> anc_bases, anc_aminos;
            Induced const sub = members.size() > 1 ? Induce(tree, member_taxids) : Induced{};
            if (sub.Leaves() >= 2) {
                size_t ignored = 0;
                anc_bases = Fitch(sub, member_bases, length, kBaseStates, ignored);
                anc_aminos = Fitch(sub, member_aminos, width, kAminoStates, ignored);
                stats.genus_trees++;
            } else {
                anc_bases = bases[members[0]];
                anc_aminos = aminos[members[0]];
            }
            stats.genus_ancestors++;
            for (size_t c = 0; c < length; c++) if (anc_bases[c] < 4) among_votes[c][anc_bases[c]]++;
            for (size_t c = 0; c < width; c++) if (anc_aminos[c] < kAminoAcids - 1) among_aa[c][anc_aminos[c]]++;
        }
        // The family's ancestral bases: Fitch over all its species.
        std::vector<uint32_t> all_taxids;
        std::vector<std::vector<uint8_t>> all_bases;
        for (size_t const i : used) {
            all_taxids.push_back(taxids[i]);
            all_bases.push_back(bases[i]);
        }
        Induced const fam_tree = Induce(tree, all_taxids);
        size_t ambiguous = 0;
        std::vector<uint8_t> const fam_anc = fam_tree.Leaves() >= 2 ? Fitch(fam_tree, all_bases, length, kBaseStates, ambiguous)
                                                                     : std::vector<uint8_t>(length, kNoBase);
        stats.family_columns += length;
        stats.family_ambiguous += ambiguous;

        Family row;
        row.family = family;
        row.gene = gene;
        row.reference = taxids[used.front()];
        row.genera = static_cast<uint16_t>(std::min<size_t>(genera.size(), UINT16_MAX));
        row.within.assign(length, kNoCode);
        row.among.assign(length, kNoCode);
        row.aa.assign(width, kNoCode);
        row.consensus.assign(length, kNoBase);
        for (size_t c = 0; c < length; c++) {
            if (within_n[c] > 0) row.within[c] = CodeOfPooled(within_sum[c] / within_n[c], within_compared[c]);
            uint32_t const n = among_votes[c][0] + among_votes[c][1] + among_votes[c][2] + among_votes[c][3];
            if (n > 0) row.among[c] = CodeOf(detail::Best(among_votes[c]).first, n);
            if (fam_anc[c] < 4 && genera_with_base[c] >= kMinConsensusGenera) {
                row.consensus[c] = fam_anc[c];
                stats.polarised++;
            }
        }
        for (size_t c = 0; c < width; c++) {
            double sum = 0;
            size_t parts = 0, compared = 0;
            if (aa_n[c] > 0) {
                sum += aa_sum[c] / aa_n[c];
                parts++;
                compared += aa_compared[c];
            }
            uint32_t n = 0;
            for (size_t k = 0; k + 1 < kAminoAcids; k++) n += among_aa[c][k];
            if (n >= kMinCompared) {
                std::array<uint32_t, kAminoAcids - 1> v{};
                for (size_t k = 0; k + 1 < kAminoAcids; k++) v[k] = among_aa[c][k];
                sum += static_cast<double>(detail::Best(v).first) / static_cast<double>(n);
                parts++;
                compared += n;
            }
            if (parts > 0) row.aa[c] = CodeOfPooled(sum / static_cast<double>(parts), compared);
        }
        uint32_t const index = static_cast<uint32_t>(families.size());
        families.push_back(std::move(row));
        for (size_t const i : used) {
            if (runs[i].empty()) continue;
            copies.emplace_back(taxids[i], gene, index, std::move(runs[i]));
        }
        return true;
    }

    // One gene's copies (as ScanGene takes them) with their rows (`rows`: taxid -> aligned row, of one width), one
    // family per task on `threads` threads; the rows in family order (deterministic).
    inline void ScanGeneFromMsa(uint32_t gene, std::vector<uint32_t> const& taxids, std::vector<std::string> const& seqs,
                                std::vector<uint32_t> const& genus_of, std::vector<uint32_t> const& family_of,
                                std::unordered_map<uint32_t, std::string> const& rows, size_t width, Tree const& tree,
                                int threads, MsaStats& stats, std::vector<Family>& families,
                                std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>>& copies) {
        threads = std::max(threads, 1);
        std::map<uint32_t, std::vector<size_t>> by_family;
        for (size_t i = 0; i < taxids.size(); i++) {
            if (family_of[i] != 0 && genus_of[i] != 0) by_family[family_of[i]].push_back(i);
        }
        std::vector<std::pair<uint32_t, std::vector<size_t>>> tasks(by_family.begin(), by_family.end());
        if (tasks.empty()) return;
        std::vector<std::vector<Family>> fam_out(tasks.size());
        std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>>> copy_out(tasks.size());
        std::vector<MsaStats> stat_out(tasks.size());
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (int64_t t = 0; t < static_cast<int64_t>(tasks.size()); t++) {
            auto const& [family, members] = tasks[static_cast<size_t>(t)];
            std::vector<uint32_t> ids, genus;
            std::vector<std::string> sequences;
            std::vector<std::string const*> member_rows;
            for (size_t const i : members) {
                ids.push_back(taxids[i]);
                genus.push_back(genus_of[i]);
                sequences.push_back(seqs[i]);
                auto const it = rows.find(taxids[i]);
                member_rows.push_back(it == rows.end() ? nullptr : &it->second);
            }
            BuildFamilyFromMsa(family, gene, ids, sequences, genus, member_rows, width, tree, stat_out[static_cast<size_t>(t)],
                               fam_out[static_cast<size_t>(t)], copy_out[static_cast<size_t>(t)]);
        }
        for (size_t t = 0; t < tasks.size(); t++) {
            uint32_t const base = static_cast<uint32_t>(families.size());
            for (auto& f : fam_out[t]) families.push_back(std::move(f));
            for (auto& [taxid, g, row, r] : copy_out[t]) copies.emplace_back(taxid, g, base + row, std::move(r));
            stats.Add(stat_out[t]);
        }
    }
}
