"""GTDB lineages for the training scripts: read from a protal database's internal_taxonomy.dmp or from
GTDB lineage strings (d__...;p__...;...;s__...), and how closely two taxa are related.

A lineage is a dict {rank: name}, names with their GTDB prefix (p__Pseudomonadota), ranks from RANKS.
"""

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
PREFIXES = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}


def from_string(lineage):
    """The lineage of a GTDB lineage string."""
    out = {}
    for token in lineage.split(";"):
        token = token.strip()
        if len(token) > 3 and token[1:3] == "__" and token[0] in PREFIXES:
            out[PREFIXES[token[0]]] = token
    return out


def from_taxonomy(path):
    """{node id: lineage} for every node of internal_taxonomy.dmp (id, parent_id, ..., name, rank, ...), and
    {name: node id}."""
    parent, name, rank = {}, {}, {}
    with open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        i_id, i_parent, i_name, i_rank = (header.index(c) for c in ("id", "parent_id", "name", "rank"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(i_id, i_parent, i_name, i_rank):
                parent[f[i_id]], name[f[i_id]], rank[f[i_id]] = f[i_parent], f[i_name], f[i_rank]
    lineages = {}

    def lineage(node):
        if node in lineages:
            return lineages[node]
        own = {rank[node]: name[node]} if rank.get(node) in RANKS else {}
        up = parent.get(node)
        result = {**(lineage(up) if up is not None and up != node and up in parent else {}), **own}
        lineages[node] = result
        return result

    for node in parent:
        lineage(node)
    return lineages, {n: i for i, n in name.items()}


def shared_rank(a, b):
    """The deepest rank at which lineages a and b agree ("species" for one species); None if not even the
    domain."""
    deepest = None
    for r in RANKS:
        if r in a and a[r] == b.get(r):
            deepest = r
        else:
            break
    return deepest
