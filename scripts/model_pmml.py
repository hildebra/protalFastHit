"""Write a scikit-learn random forest as the PMML model protal scores taxa with, and score such a
model in Python the way protal's cPMML does.

protal reads the model with cPMML (C++). sklearn2pmml, which wrote the models before, runs the Java
JPMML converter; this module writes the same kind of model directly, so training needs no Java: a
MiningModel whose Segmentation averages one TreeModel per tree, as predict_proba does.

protal's probabilities equal sklearn's bit for bit:
- sklearn compares a feature as float32 with a double threshold, cPMML compares doubles. Each
  threshold is written as the largest double x for which float32(x) <= threshold, so both send
  every value down the same branch (FLOAT32_SPLIT).
- A leaf's probability in cPMML is its TRUE count over the sum of its counts. The counts written are
  1 - p and p, which sum to exactly 1.0 in doubles, so cPMML gets sklearn's p unchanged.
- cPMML sums the trees in file order and divides by their number, as predict_proba does with
  n_jobs=1 (with more jobs, the summation order and so the last bit may differ).
"""

import xml.etree.ElementTree as ET
from xml.sax.saxutils import quoteattr

import numpy as np

LABELS = ("FALSE", "TRUE")

# protal warns when it loads a model with this annotation (profiler::kPlaceholderModelMarker).
PLACEHOLDER_MARKER = "protal:placeholder"
# Database member of each read type's model (kReadTypes in src/ReadType.h).
MODEL_FILES = {"pe": "model_pe.xml", "se": "model_se.xml", "pb": "model_PB.xml", "ont": "model_ONT.xml"}


def write_placeholder(path, read_type):
    """An untrained model for `read_type`: it scores every taxon 0, so protal reports no species (with
    --knob 0, every taxon with reads). protal warns when it loads it."""
    field = "fragments"
    schema = (f'<MiningSchema><MiningField name="truth" usageType="predicted"/><MiningField name="{field}"/>'
              '</MiningSchema>')
    text = (
        '<?xml version="1.0" encoding="UTF-8"?>\n<PMML version="4.4">\n'
        f' <Header description="protal placeholder model for read type {read_type}: untrained, scores every taxon 0">\n'
        '  <Application name="protal scripts/placeholder_models.py"/>\n'
        f'  <Annotation>{PLACEHOLDER_MARKER}</Annotation>\n'
        f'  <Annotation>replace with a trained model: protal --add_model MODEL --read_type {read_type} --db DB</Annotation>\n'
        ' </Header>\n'
        ' <DataDictionary numberOfFields="2">\n'
        f'  <DataField name="truth" optype="categorical" dataType="string"><Value value="{LABELS[0]}"/>'
        f'<Value value="{LABELS[1]}"/></DataField>\n'
        f'  <DataField name="{field}" optype="continuous" dataType="double"/>\n'
        ' </DataDictionary>\n'
        f' <MiningModel functionName="classification">\n  {schema}\n  <Segmentation multipleModelMethod="average">\n'
        f'   <Segment id="1"><True/><TreeModel functionName="classification" splitCharacteristic="binarySplit">{schema}\n'
        f'<Node score="{LABELS[0]}" recordCount="1.0"><True/><ScoreDistribution value="{LABELS[0]}" recordCount="1.0"/>'
        f'<ScoreDistribution value="{LABELS[1]}" recordCount="0.0"/></Node>'
        '</TreeModel></Segment>\n  </Segmentation>\n </MiningModel>\n</PMML>\n')
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def float32_split(threshold):
    """The largest double x with float32(x) <= threshold: `x <= float32_split(t)` is sklearn's
    `float32(x) <= t` for every double x."""
    t = float(threshold)
    below = np.float32(t)
    if float(below) > t:
        below = np.nextafter(below, np.float32(-np.inf))
    above = np.nextafter(below, np.float32(np.inf))
    if not np.isfinite(above):
        return t
    middle = (float(below) + float(above)) / 2.0  # exact: a float32 midpoint fits in a double
    # A double exactly halfway rounds to the float32 with the even significand.
    even = int(np.array(below, dtype=np.float32).view(np.uint32)) & 1 == 0
    return middle if even else float(np.nextafter(middle, -np.inf))


def write_forest(forest, features, path, annotations=()):
    """Write a fitted RandomForestClassifier with classes 0/1 (or False/True), trained on the
    columns `features` in this order, as PMML for protal. `annotations` go into the header."""
    try:
        classes = [int(c) for c in forest.classes_]
    except (TypeError, ValueError):
        classes = None
    if classes != [0, 1]:
        raise ValueError(f"expected a binary forest with classes 0 and 1, got {list(forest.classes_)}")
    if forest.n_features_in_ != len(features):
        raise ValueError(f"the forest has {forest.n_features_in_} inputs, not {len(features)} features")
    out = []
    w = out.append
    w('<?xml version="1.0" encoding="UTF-8"?>\n<PMML version="4.4">\n')
    w(' <Header description="protal presence model: probability that a taxon is present">\n')
    w('  <Application name="protal scripts/random_forest_cmdline.py"/>\n')
    for note in annotations:
        w(f'  <Annotation>{_text(note)}</Annotation>\n')
    w(' </Header>\n')
    w(f' <DataDictionary numberOfFields="{len(features) + 1}">\n')
    w(f'  <DataField name="truth" optype="categorical" dataType="string">'
      f'<Value value="{LABELS[0]}"/><Value value="{LABELS[1]}"/></DataField>\n')
    for f in features:
        w(f'  <DataField name={quoteattr(f)} optype="continuous" dataType="double"/>\n')
    w(' </DataDictionary>\n')
    schema = ('<MiningSchema><MiningField name="truth" usageType="predicted"/>'
              + "".join(f'<MiningField name={quoteattr(f)}/>' for f in features) + '</MiningSchema>')
    w(f' <MiningModel functionName="classification">\n  {schema}\n  <Segmentation multipleModelMethod="average">\n')
    names = [quoteattr(f) for f in features]
    for i, estimator in enumerate(forest.estimators_):
        tree = estimator.tree_
        left, right = tree.children_left, tree.children_right
        feature, threshold = tree.feature, tree.threshold
        count = tree.weighted_n_node_samples
        probability = leaf_probabilities(tree)
        w(f'   <Segment id="{i + 1}"><True/>'
          f'<TreeModel functionName="classification" splitCharacteristic="binarySplit">{schema}\n')
        # Depth first, left child first, without recursion (trees without a leaf limit get deep).
        stack = [(0, "<True/>", False)]
        while stack:
            node, predicate, closing = stack.pop()
            if closing:
                w('</Node>')
                continue
            p = float(probability[node])
            w(f'<Node score="{LABELS[p > 0.5]}" recordCount="{float(count[node])!r}">{predicate}')
            if left[node] == -1:
                w(f'<ScoreDistribution value="{LABELS[0]}" recordCount="{1.0 - p!r}"/>'
                  f'<ScoreDistribution value="{LABELS[1]}" recordCount="{p!r}"/></Node>')
                continue
            split = float32_split(threshold[node])
            name = names[feature[node]]
            stack.append((node, None, True))
            stack.append((right[node], f'<SimplePredicate field={name} operator="greaterThan" value="{split!r}"/>', False))
            stack.append((left[node], f'<SimplePredicate field={name} operator="lessOrEqual" value="{split!r}"/>', False))
        w('</TreeModel></Segment>\n')
    w('  </Segmentation>\n </MiningModel>\n</PMML>\n')
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("".join(out))


def leaf_probabilities(tree):
    """P(TRUE) of every node as the tree's predict_proba gives it. scikit-learn 1.4 and later store
    class fractions in tree_.value and return them as they are; older versions store weighted
    counts and divide by their sum."""
    value = tree.value[:, 0, :]
    total = value.sum(axis=1)
    if np.allclose(total, 1.0, rtol=0, atol=1e-9):
        return value[:, 1].copy()
    return value[:, 1] / np.where(total == 0, 1.0, total)


def _text(value):
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class PmmlForest:
    """A PMML forest as write_forest writes it, scored as cPMML scores it: the probability of TRUE
    is the average over the trees of each leaf's TRUE count over the sum of its counts."""

    def __init__(self, path):
        root = ET.parse(path).getroot()
        for elem in root.iter():
            if "}" in elem.tag:
                elem.tag = elem.tag.split("}", 1)[1]
        model = root.find("MiningModel")
        if model is None:
            raise ValueError(f"{path}: no MiningModel")
        schema = model.find("MiningSchema")
        self.features = [f.get("name") for f in schema.findall("MiningField")
                         if f.get("usageType", "active") == "active"]
        index = {name: i for i, name in enumerate(self.features)}
        segmentation = model.find("Segmentation")
        if segmentation.get("multipleModelMethod") != "average":
            raise ValueError(f"{path}: the trees are combined by {segmentation.get('multipleModelMethod')}, not averaged")
        self.trees = [self._compile(segment.find("TreeModel").find("Node"), index)
                      for segment in segmentation.findall("Segment")]

    @staticmethod
    def _compile(top, index):
        feature, threshold, left, right, leaf = [], [], [], [], []

        def add():
            feature.append(-1)
            threshold.append(np.nan)
            left.append(-1)
            right.append(-1)
            leaf.append(np.nan)
            return len(feature) - 1

        stack = [(top, add())]
        while stack:
            node, n = stack.pop()
            children = node.findall("Node")
            if not children:
                counts = {d.get("value"): float(d.get("recordCount")) for d in node.findall("ScoreDistribution")}
                total = 0.0
                for d in node.findall("ScoreDistribution"):
                    total += float(d.get("recordCount"))
                leaf[n] = counts.get("TRUE", 0.0) / total
                continue
            if len(children) != 2:
                raise ValueError("expected binary splits")
            first, second = (c.find("SimplePredicate") for c in children)
            if (first.get("operator"), second.get("operator")) != ("lessOrEqual", "greaterThan") or \
                    first.get("field") != second.get("field") or first.get("value") != second.get("value"):
                raise ValueError("expected splits 'x <= t' then 'x > t' on one field")
            feature[n] = index[first.get("field")]
            threshold[n] = float(first.get("value"))
            left[n], right[n] = add(), add()
            stack += [(children[1], right[n]), (children[0], left[n])]
        return (np.array(feature), np.array(threshold), np.array(left), np.array(right), np.array(leaf))

    def predict(self, X):
        """Probability of TRUE for the rows of X (a DataFrame with the model's features, or an
        array with them in model order), as doubles, compared as doubles."""
        if hasattr(X, "columns"):
            X = X[self.features].to_numpy(dtype=np.float64)
        X = np.asarray(X, dtype=np.float64)
        rows = np.arange(len(X))
        total = np.zeros(len(X))
        for feature, threshold, left, right, leaf in self.trees:
            node = np.zeros(len(X), dtype=np.int64)
            active = feature[node] >= 0
            while active.any():
                at = node[active]
                go_left = X[rows[active], feature[at]] <= threshold[at]
                node[active] = np.where(go_left, left[at], right[at])
                active = feature[node] >= 0
            total += leaf[node]
        return total / len(self.trees)

    def node_count(self):
        return sum(len(tree[0]) for tree in self.trees)
