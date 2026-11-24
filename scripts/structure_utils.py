"""Small auditable Q-matching and Newick helpers, independent of plotting."""
from dataclasses import dataclass, field
import re
import gzip
import numpy as np


def best_assignment(score):
    """Maximize one-to-one rectangular assignment for rows <= columns.

    Exact subset dynamic programming, sufficient for the tested K <= 8.
    Returns one selected column for each row.
    """
    score = np.asarray(score, dtype=float)
    n, m = score.shape
    assert n <= m
    states = {0: (0., ())}
    for row in range(n):
        updated = {}
        for mask, (value, assignment) in states.items():
            for col in range(m):
                if mask & (1 << col):
                    continue
                new_mask = mask | (1 << col)
                candidate = (value + score[row, col], assignment + (col,))
                if new_mask not in updated or candidate[0] > updated[new_mask][0]:
                    updated[new_mask] = candidate
        states = updated
    return list(max(states.values(), key=lambda x: x[0])[1])


def align_q(reference, candidate):
    perm = best_assignment(reference.T @ candidate)
    matched = candidate[:, perm]
    return matched, perm, float(np.sqrt(np.mean((reference-matched)**2)))


def read_result(path, expected_ids):
    text = path.read_text(errors="replace")
    number = r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
    k = int(re.search(r"(\d+) populations assumed", text).group(1))
    values = {}
    for name, pattern in {
        "ln_probability": rf"Estimated Ln Prob of Data\s*=\s*{number}",
        "mean_ln_likelihood": rf"Mean value of ln likelihood\s*=\s*{number}",
        "variance_ln_likelihood": rf"Variance of ln likelihood\s*=\s*{number}",
        "alpha": rf"Mean value of alpha\s*=\s*{number}",
    }.items():
        match = re.search(pattern, text)
        values[name] = float(match.group(1)) if match else float("nan")
    qrows = {}
    section = text.split("Inferred ancestry of individuals:", 1)[1]
    for line in section.splitlines():
        match = re.match(r"\s*\d+\s+(\S+)\s+\([^)]*\)\s*:\s*(.*)", line)
        if match:
            row = [float(value) for value in match.group(2).split()[:k]]
            if len(row) == k:
                qrows[match.group(1)] = row
    assert set(qrows) == set(expected_ids), (path, set(expected_ids)-set(qrows))
    raw = np.array([qrows[name] for name in expected_ids])
    assert np.all((raw >= 0) & (raw <= 1))
    assert np.max(np.abs(raw.sum(axis=1)-1)) <= .005
    # Correct only the unavoidable 3-decimal output rounding.
    q = raw / raw.sum(axis=1, keepdims=True)
    return k, values, q


def read_likelihood_trace(path, burnin, numreps):
    values = []
    text = (gzip.open(path,"rt",errors="replace").read() if path.suffix == ".gz" else path.read_text(errors="replace")).split("Finished initialization; starting MCMC", 1)[1]
    for line in text.splitlines():
        if re.fullmatch(r"\s*[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?\s*", line):
            values.append(float(line))
    expected = burnin+numreps
    # At K=1 the two printed membership summaries each add a standalone 1.000.
    assert expected <= len(values) <= expected+2, (path, len(values), expected)
    assert all(value < 0 for value in values[:expected])
    assert all(value >= 0 for value in values[expected:])
    return np.array(values[burnin:expected])


def split_rhat(traces):
    n = min(map(len, traces)) // 2
    chains = np.array([half for t in traces for half in [t[:n], t[-n:]]])
    within = np.mean(np.var(chains, axis=1, ddof=1))
    between = n * np.var(chains.mean(axis=1), ddof=1)
    return float(np.sqrt(((n-1)/n*within + between/n)/within)) if within > 0 else float("nan")


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.
    children: list = field(default_factory=list)


def read_newick(path):
    tokens = re.findall(r"\(|\)|,|:|;|[^\s(),:;]+", path.read_text())
    cursor = 0
    def parse():
        nonlocal cursor
        node = Node()
        if tokens[cursor] == "(":
            cursor += 1
            while True:
                node.children.append(parse())
                if tokens[cursor] == ",":
                    cursor += 1
                else:
                    assert tokens[cursor] == ")"
                    cursor += 1
                    break
        if cursor < len(tokens) and tokens[cursor] not in (":", ",", ")", ";"):
            node.name = tokens[cursor].strip("'")
            cursor += 1
        if cursor < len(tokens) and tokens[cursor] == ":":
            cursor += 1
            node.length = float(tokens[cursor])
            cursor += 1
        return node
    root = parse()
    assert tokens[cursor] == ";"
    return root


def leaves(node):
    return [node] if not node.children else [leaf for child in node.children for leaf in leaves(child)]


def normalize_tree_labels(root):
    mapping = []
    for tip in leaves(root):
        original = tip.name
        match = re.search(r"(F\d{5})$", original)
        if match:
            tip.name = match.group(1)
        mapping.append((original, tip.name))
    assert len({tip.name for tip in leaves(root)}) == len(mapping)
    return mapping


def display_root_on_tip(root, name):
    """Reorient the unrooted graph at a tip solely for rectangular display."""
    adjacency = {}
    def visit(node, parent=None):
        adjacency.setdefault(node, [])
        if parent is not None:
            adjacency[node].append((parent, node.length))
            adjacency[parent].append((node, node.length))
        for child in node.children:
            visit(child, node)
    visit(root)
    tip = next(n for n in adjacency if not n.children and n.name == name)
    neighbor, length = adjacency[tip][0]
    def copy(node, previous, edge_length):
        children = [copy(n, node, l) for n, l in adjacency[node] if n is not previous]
        label = node.name if not children else ""
        return Node(label, edge_length, children)
    return Node("", 0., [copy(tip, neighbor, length/2), copy(neighbor, tip, length/2)])
