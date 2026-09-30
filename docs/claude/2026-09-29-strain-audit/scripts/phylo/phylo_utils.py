"""Small, dependency-free tree utilities (Newick I/O, splits, RF, patristic distances).

Trees are nested Node objects. Splits are frozensets of leaf names on the side that does
not contain a fixed anchor leaf, so they are comparable between unrooted trees on the same
leaf set.
"""
import itertools
import math
import re


class Node:
    __slots__ = ("name", "length", "children", "parent", "label")

    def __init__(self, name=None, length=None):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None
        self.label = None  # internal node label (e.g. bootstrap support)

    def add(self, child):
        child.parent = self
        self.children.append(child)
        return child

    def is_leaf(self):
        return not self.children

    def leaves(self):
        if self.is_leaf():
            return [self]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out

    def nodes(self):
        out = [self]
        for c in self.children:
            out.extend(c.nodes())
        return out


def parse_newick(text):
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    pos = 0
    n = len(text)

    def parse_name():
        nonlocal pos
        if pos < n and text[pos] == "'":
            end = text.index("'", pos + 1)
            name = text[pos + 1:end]
            pos = end + 1
            return name
        start = pos
        while pos < n and text[pos] not in ",():;[":
            pos += 1
        return text[start:pos].strip()

    def parse_len():
        nonlocal pos
        if pos < n and text[pos] == ":":
            pos += 1
            start = pos
            while pos < n and text[pos] not in ",();[":
                pos += 1
            return float(text[start:pos])
        return None

    def skip_comment():
        nonlocal pos
        if pos < n and text[pos] == "[":
            pos = text.index("]", pos) + 1

    def parse_subtree():
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                node.add(parse_subtree())
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"bad newick at {pos}: {text[pos:pos+20]}")
            lab = parse_name()
            node.label = lab if lab != "" else None
        else:
            node.name = parse_name()
        skip_comment()
        node.length = parse_len()
        skip_comment()
        return node

    root = parse_subtree()
    return root


def to_newick(node, digits=8):
    def rec(nd):
        if nd.is_leaf():
            s = nd.name
        else:
            s = "(" + ",".join(rec(c) for c in nd.children) + ")"
            if nd.label is not None:
                s += str(nd.label)
        if nd.length is not None:
            s += ":" + f"{nd.length:.{digits}g}"
        return s
    return rec(node) + ";"


def read_tree(path):
    with open(path) as fh:
        return parse_newick(fh.read())


def leaf_names(root):
    return sorted(l.name for l in root.leaves())


def prune(root, keep):
    """Return a copy of the tree restricted to leaves in `keep`, suppressing unary nodes."""
    keep = set(keep)

    def rec(nd):
        if nd.is_leaf():
            if nd.name in keep:
                c = Node(nd.name, nd.length)
                return c
            return None
        kids = [k for k in (rec(c) for c in nd.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            k = kids[0]
            if nd.length is not None or k.length is not None:
                k.length = (k.length or 0.0) + (nd.length or 0.0)
            return k
        c = Node(None, nd.length)
        c.label = nd.label
        for k in kids:
            c.add(k)
        return c

    r = rec(root)
    r.length = None
    # an unrooted tree read as a rooted binary tree may leave a degree-2 root: fine for splits
    return r


def splits(root, leafset=None, with_info=False):
    """Non-trivial splits of the (unrooted) tree. Each split is the frozenset of leaves on the
    side NOT containing the lexicographically smallest leaf. Returns a dict
    split -> (branch_length, support) if with_info else a set."""
    all_leaves = frozenset(l.name for l in root.leaves())
    if leafset is not None:
        assert frozenset(leafset) == all_leaves, (sorted(leafset), sorted(all_leaves))
    anchor = min(all_leaves)
    n = len(all_leaves)
    out = {}

    def rec(nd):
        if nd.is_leaf():
            return frozenset([nd.name])
        s = frozenset()
        for c in nd.children:
            s = s | rec(c)
        if nd is not root:
            side = s if anchor not in s else all_leaves - s
            if 1 < len(side) < n - 1:
                sup = None
                if nd.label is not None:
                    try:
                        sup = float(str(nd.label).split("/")[-1])
                    except ValueError:
                        sup = None
                bl = nd.length or 0.0
                if side in out:  # degree-2 root: both root edges give the same split
                    obl, osup = out[side]
                    out[side] = (obl + bl, osup if osup is not None else sup)
                else:
                    out[side] = (bl, sup)
        return s

    rec(root)
    return out if with_info else set(out)


def rf(t1, t2):
    """Robinson-Foulds distance on the common leaf set; returns (rf, max_rf, n_splits_t1, n_splits_t2)."""
    common = set(leaf_names(t1)) & set(leaf_names(t2))
    a = prune(t1, common)
    b = prune(t2, common)
    s1 = splits(a)
    s2 = splits(b)
    d = len(s1 ^ s2)
    n = len(common)
    max_rf = 2 * (n - 3) if n > 3 else 0
    return d, max_rf, len(s1), len(s2)


def patristic(root):
    """Dict (a,b) -> path length between leaves."""
    # distances from each node to its leaves
    leaves = root.leaves()
    # compute root-to-node depth, then use LCA via parent chains
    depth = {}

    def rec(nd, d):
        depth[id(nd)] = d
        for c in nd.children:
            rec(c, d + (c.length or 0.0))

    rec(root, 0.0)
    anc = {}
    for l in leaves:
        chain = []
        x = l
        while x is not None:
            chain.append(x)
            x = x.parent
        anc[l.name] = chain
    out = {}
    for a, b in itertools.combinations(leaves, 2):
        sa = {id(x) for x in anc[a.name]}
        lca = next(x for x in anc[b.name] if id(x) in sa)
        d = depth[id(a)] + depth[id(b)] - 2 * depth[id(lca)]
        out[(a.name, b.name)] = d
        out[(b.name, a.name)] = d
    return out


def unroot(root):
    """Turn a degree-2 root into a trifurcation (so terminal edges are unrooted edges)."""
    if len(root.children) == 2:
        a, b = root.children
        if a.is_leaf() and b.is_leaf():
            return root
        if a.is_leaf():
            a, b = b, a
        # a is internal: remove it, give its length to b
        b.length = (b.length or 0.0) + (a.length or 0.0)
        root.children = [b] + a.children
        for c in root.children:
            c.parent = root
        if a.label is not None and root.label is None:
            root.label = None
    return root


def fix_parents(root):
    for nd in root.nodes():
        for c in nd.children:
            c.parent = nd
    root.parent = None
    return root


def terminal_lengths(root):
    return {l.name: (l.length or 0.0) for l in root.leaves()}


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2.0
            i = j + 1
        return r
    return pearson(ranks(x), ranks(y))


def quartet_distance(t1, t2, max_quartets=None):
    """Fraction of 4-leaf subsets whose induced topology differs (resolved vs resolved,
    unresolved counted as differing). Exact enumeration (fine for <= ~20 leaves)."""
    common = sorted(set(leaf_names(t1)) & set(leaf_names(t2)))
    s1 = splits(prune(t1, common))
    s2 = splits(prune(t2, common))
    allset = frozenset(common)

    def qtop(sp, q):
        qs = frozenset(q)
        for s in sp:
            inter = qs & s
            if len(inter) == 2:
                return frozenset([frozenset(inter), frozenset(qs - inter)])
        return None

    diff = 0
    tot = 0
    for q in itertools.combinations(common, 4):
        tot += 1
        a = qtop(s1, q)
        b = qtop(s2, q)
        if a is None or b is None or a != b:
            diff += 1
    return diff / tot if tot else float("nan")
