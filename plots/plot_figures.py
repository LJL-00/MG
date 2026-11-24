#!/usr/bin/env python3
"""Draw manuscript panels from user-supplied analysis outputs.

Inputs are CSV/TSV matrices with sample or gene IDs in the first column.
The integrated figure uses an externally supplied Newick tree and annotation
tables. No tree inference, clade assignment or phenotype classification is
performed here. The RNA panel uses supplied differential-expression results.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Rectangle


CLINICAL = "#C44E52"
NONCLINICAL = "#4C72B0"
MISSING = "#D9D9D9"


def configure():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
        "axes.labelsize": 8, "axes.titlesize": 9, "axes.linewidth": .7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none"})


def table(path, *, numeric=False):
    frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    frame.index = frame.index.map(str)
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError(f"Duplicate identifiers in {path}")
    if numeric:
        frame = frame.apply(pd.to_numeric, errors="raise")
        if np.isinf(frame.to_numpy(float)).any():
            raise ValueError(f"Infinite values in {path}")
    return frame


def ids_file(path):
    result = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()
              if line.strip() and not line.lstrip().startswith("#")]
    if len(result) != len(set(result)):
        raise ValueError("Duplicate ID in selection list")
    return result


def export(fig, path, dpi=600):
    path = Path(path)
    if path.suffix.lower() not in {".pdf", ".svg", ".png", ".tif", ".tiff"}:
        raise ValueError("Output extension must be pdf, svg, png, tif or tiff")
    path.parent.mkdir(parents=True, exist_ok=True)
    extra = {"pil_kwargs": {"compression": "tiff_lzw"}} if path.suffix.lower() in {".tif", ".tiff"} else {}
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white", **extra)
    plt.close(fig)
    return path


def clean_axis(ax):
    for edge in ["top", "right"]:
        ax.spines[edge].set_visible(False)


def source_color(value):
    name = str(value).strip().lower()
    return CLINICAL if name == "clinical" else NONCLINICAL if name in {
        "environment", "environmental", "nonclinical", "non-clinical"} else "#777777"


def square_matrix(frame, *, symmetry=False):
    if set(frame.index) != set(frame.columns):
        raise ValueError("A square matrix must have the same row and column IDs")
    frame = frame.loc[frame.index, frame.index]
    if symmetry and not np.allclose(frame, frame.T, equal_nan=True, atol=1e-7):
        raise ValueError("Expected a symmetric matrix")
    return frame


def matrix_panel(ax, frame, title, *, vmin, vmax, cmap="viridis", annotate=False):
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad(MISSING)
    im = ax.imshow(np.ma.masked_invalid(frame.to_numpy(float)), vmin=vmin,
                   vmax=vmax, cmap=palette, interpolation="nearest", aspect="equal")
    ax.set_xticks(range(len(frame.columns)), frame.columns, rotation=90)
    ax.set_yticks(range(len(frame.index)), frame.index)
    ax.tick_params(length=0, pad=2)
    ax.set_title(title, loc="left", fontweight="bold")
    if annotate and max(frame.shape) <= 20:
        for y in range(len(frame)):
            for x in range(len(frame.columns)):
                value = frame.iloc[y, x]
                if np.isfinite(value):
                    ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=6)
    return im


def ani_figure(matrix, subset=None, *, vmin=86., vmax=100.):
    matrix = square_matrix(matrix, symmetry=True)
    finite = matrix.to_numpy(float)[np.isfinite(matrix.to_numpy(float))]
    if not len(finite) or (finite < 0).any() or (finite > 100).any():
        raise ValueError("ANI values must be percentages between 0 and 100")
    if not vmin < vmax:
        raise ValueError("ANI colour limits must increase")
    panels = [matrix]
    if subset:
        if not set(subset).issubset(matrix.index):
            raise ValueError("ANI subset contains unknown samples")
        panels.append(matrix.loc[subset, subset])
    fig, axes = plt.subplots(1, len(panels), figsize=(7.2 if len(panels) == 2 else 5.2, 5),
                             squeeze=False, layout="constrained")
    for i, frame in enumerate(panels):
        im = matrix_panel(axes[0, i], frame, "a" if i == 0 else "b", vmin=vmin,
                          vmax=vmax, cmap="viridis", annotate=False)
    fig.colorbar(im, ax=list(axes[0]), fraction=.035, pad=.025, label="Reciprocal mean ANI (%)")
    return fig


def pca_panel(ax, coordinates, metadata, variance=None, *, flip_pc1=False):
    if not {"PC1", "PC2"}.issubset(coordinates.columns):
        raise ValueError("PCA table needs PC1 and PC2")
    if not set(coordinates.index).issubset(metadata.index) or "group" not in metadata:
        raise ValueError("Metadata need a group column and all PCA sample IDs")
    x = coordinates["PC1"].to_numpy(float) * (-1 if flip_pc1 else 1)
    y = coordinates["PC2"].to_numpy(float)
    if not np.isfinite(np.c_[x, y]).all():
        raise ValueError("PCA coordinates must be finite")
    groups = metadata.loc[coordinates.index, "group"].astype(str)
    for group in dict.fromkeys(groups):
        mask = groups.to_numpy() == group
        ax.scatter(x[mask], y[mask], s=32, color=source_color(group), edgecolor="white",
                   linewidth=.4, label=group, zorder=3)
    for i, name in enumerate(coordinates.index):
        ax.annotate(name, (x[i], y[i]), xytext=(4, 3 if i % 2 else -9),
                    textcoords="offset points", fontsize=6)
    suffix = ["", ""]
    if variance is not None:
        lookup = dict(zip(variance["PC"].astype(int), variance["variance_fraction"].astype(float)))
        suffix = [f" ({lookup[i]:.2%})" for i in [1, 2]]
    ax.set_xlabel("PC1" + suffix[0]); ax.set_ylabel("PC2" + suffix[1])
    ax.legend(frameon=False, fontsize=6, loc="best")
    ax.margins(.15)
    clean_axis(ax)


def differential_categories(frame, *, effect="log2FoldChange", adjusted_p="padj", alpha=.05, fold_change=1.):
    values = frame[[effect, adjusted_p]].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(values).all().all() or not values[adjusted_p].between(0, 1).all():
        raise ValueError("Differential-expression effects must be finite and adjusted P values in [0,1]")
    if not 0 < alpha < 1 or fold_change <= 0:
        raise ValueError("Invalid differential-expression display threshold")
    labels = np.full(len(values), "Other", dtype=object)
    labels[(values[adjusted_p] < alpha) & (values[effect] >= fold_change)] = "Clinical higher"
    labels[(values[adjusted_p] < alpha) & (values[effect] <= -fold_change)] = "Nonclinical higher"
    return values, labels


def volcano_panel(ax, frame, *, effect="log2FoldChange", adjusted_p="padj", alpha=.05, fold_change=1.):
    values, labels = differential_categories(frame, effect=effect, adjusted_p=adjusted_p,
                                              alpha=alpha, fold_change=fold_change)
    positive = values.loc[values[adjusted_p] > 0, adjusted_p]
    floor = max(float(positive.min()) / 10, np.finfo(float).tiny) if len(positive) else np.finfo(float).tiny
    height = -np.log10(values[adjusted_p].clip(lower=floor))
    for name, color in [("Other", "#BDBDBD"), ("Nonclinical higher", NONCLINICAL), ("Clinical higher", CLINICAL)]:
        mask = labels == name
        ax.scatter(values.loc[mask, effect], height.loc[mask], s=5, color=color, linewidths=0,
                   rasterized=True, label=f"{name} ({int(mask.sum())})")
    ax.axhline(-np.log10(alpha), color="#555555", ls="--", lw=.6)
    for value in [-fold_change, fold_change]:
        ax.axvline(value, color="#555555", ls="--", lw=.6)
    ax.set_xlabel("log2 fold change (clinical / nonclinical)")
    ax.set_ylabel("-log10 adjusted P")
    ax.legend(frameon=False, fontsize=6, loc="best")
    clean_axis(ax)


def gene_set_panels(ax_nes, ax_p, results):
    needed = {"gene_set", "GSEA_NES", "GSEA_phenotype_permutation_p", "GSEA_fdr", "n_label_permutations"}
    if not needed.issubset(results):
        raise ValueError("Gene-set table must use the columns from 06_gene_sets.py")
    if results["gene_set"].duplicated().any():
        raise ValueError("Duplicate gene-set name")
    numeric = results[list(needed - {"gene_set"})].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric).all().all():
        raise ValueError("Gene-set results must be finite")
    p = numeric["GSEA_phenotype_permutation_p"]
    if not p.between(0, 1, inclusive="right").all():
        raise ValueError("Exact permutation P values must lie in (0,1]")
    y = np.arange(len(results))
    color = [CLINICAL if x > 0 else NONCLINICAL for x in numeric["GSEA_NES"]]
    ax_nes.barh(y, numeric["GSEA_NES"], color=color, height=.65)
    ax_nes.axvline(0, color="#555555", lw=.6)
    ax_nes.set_yticks(y, results.gene_set.str.replace("_", " "))
    ax_nes.invert_yaxis(); ax_nes.set_xlabel("Normalized enrichment score")
    ax_p.scatter(p, y, s=18, color=["black" if q < .05 else "#888888" for q in numeric["GSEA_fdr"]])
    ax_p.set_yticks(y, []); ax_p.set_ylim(ax_nes.get_ylim())
    counts = numeric["n_label_permutations"].unique()
    if len(counts) == 1 and counts[0] >= 2:
        for value, pattern in [(1/counts[0], ":"), (2/counts[0], "--")]:
            ax_p.axvline(value, ls=pattern, lw=.6, color="#666666")
    ax_p.set_xscale("log"); ax_p.set_xlim(min(float(p.min())*.8, .01), 1.2)
    ax_p.set_xlabel("Two-sided exact P")
    clean_axis(ax_nes); clean_axis(ax_p)


def rna_figure(coordinates, metadata, differential, results, variance=None, *, flip_pc1=False,
               effect="log2FoldChange", adjusted_p="padj"):
    fig = plt.figure(figsize=(10.4, 7), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.25])
    a, b = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    pca_panel(a, coordinates, metadata, variance, flip_pc1=flip_pc1)
    volcano_panel(b, differential, effect=effect, adjusted_p=adjusted_p)
    lower = grid[1, :].subgridspec(1, 2, width_ratios=[3, 1.3])
    c, d = fig.add_subplot(lower[0]), fig.add_subplot(lower[1])
    gene_set_panels(c, d, results)
    for ax, title in [(a, "a"), (b, "b"), (c, "c")]:
        ax.set_title(title, loc="left", fontweight="bold")
    return fig


def heatmap_figure(matrix, *, genes=None, log2=False, row_zscore=False, title="Expression"):
    if genes is not None:
        if not set(genes).issubset(matrix.index):
            raise ValueError("Gene selection contains unknown IDs")
        matrix = matrix.loc[genes]
    values = matrix.to_numpy(float)
    if not np.isfinite(values).all() or not values.size:
        raise ValueError("Expression matrix must contain finite values")
    if log2:
        if (values < 0).any():
            raise ValueError("log2(x+1) requires nonnegative input")
        values = np.log2(values + 1)
    if row_zscore:
        sd = values.std(axis=1, keepdims=True)
        values = np.divide(values-values.mean(axis=1, keepdims=True), sd,
                           out=np.zeros_like(values), where=sd > 0)
    transformed = pd.DataFrame(values, index=matrix.index, columns=matrix.columns)
    fig, ax = plt.subplots(figsize=(max(4.2, matrix.shape[1]*.45), max(2.4, len(matrix)*.18)), layout="constrained")
    limit = float(np.abs(values).max()) or 1.
    vmin, vmax = (-limit, limit) if row_zscore else (float(values.min()), float(values.max()) or 1.)
    im = matrix_panel(ax, transformed, title, vmin=vmin, vmax=vmax,
                       cmap="RdBu_r" if row_zscore else "viridis")
    ax.set_aspect("auto")
    label = "Gene-wise z score" if row_zscore else "log2(expression + 1)" if log2 else "Supplied expression value"
    fig.colorbar(im, ax=ax, fraction=.035, pad=.025, label=label)
    return fig


def variants_figure(frame, *, contig="CHROM", position="POS", impact="impact", gene="gene_id", label_genes=False):
    needed = {contig, position, impact, gene}
    if frame.empty or not needed.issubset(frame):
        raise ValueError("Candidate table needs contig, position, impact and gene-ID columns")
    if frame[contig].isna().any():
        raise ValueError("Candidate contig identifiers must be present")
    positions = pd.to_numeric(frame[position], errors="raise").to_numpy(float)
    if not np.isfinite(positions).all() or (positions < 1).any() or not np.equal(positions, np.floor(positions)).all():
        raise ValueError("Candidate positions must be positive, 1-based integers")
    contigs = frame[contig].astype(str)
    effects = frame[impact].fillna("Unknown").astype(str)
    order = list(dict.fromkeys(contigs))
    known = [name for name in ["HIGH", "MODERATE", "LOW", "MODIFIER"] if name in set(effects)]
    levels = known + sorted(set(effects) - set(known))
    fig, axes = plt.subplots(1, 2, figsize=(9, max(2.6, .32*len(order)+1.3)),
                             gridspec_kw={"width_ratios": [3, 1.5]}, layout="constrained")
    left = np.zeros(len(order))
    row_index = contigs.map({name:i for i,name in enumerate(order)}).to_numpy(int)
    for i, effect in enumerate(levels):
        selected = effects.to_numpy() == effect
        color = plt.get_cmap("tab10")(i % 10)
        axes[0].scatter(positions[selected]/1e6, row_index[selected], s=16,
                        color=color, label=effect, linewidths=0)
        counts = np.array([int((selected & (contigs.to_numpy() == name)).sum()) for name in order])
        axes[1].barh(range(len(order)), counts, left=left, height=.65, color=color)
        left += counts
    if label_genes:
        for i, value in enumerate(frame[gene]):
            if pd.notna(value):
                axes[0].annotate(str(value), (positions[i]/1e6, row_index[i]),
                                 xytext=(2, 3), textcoords="offset points", fontsize=5)
    for ax in axes:
        ax.set_ylim(len(order)-.5, -.5)
        clean_axis(ax)
    axes[0].set_yticks(range(len(order)), order)
    axes[1].set_yticks(range(len(order)), [])
    axes[0].set_xlabel("Reference position (Mb)")
    axes[1].set_xlabel("Candidate records")
    axes[0].set_title("a  Candidate positions", loc="left", fontweight="bold")
    axes[1].set_title("b  Recorded impacts", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=6)
    return fig


def qc_values(frame, *, sample="unit", isolate="isolate", mapping="mapped_fraction",
              assigned="assigned_fragments", mapping_total=None):
    needed = {sample, isolate, mapping, assigned} | ({mapping_total} if mapping_total else set())
    if frame.empty or not needed.issubset(frame):
        raise ValueError("QC table is missing a requested column")
    if frame[[sample, isolate]].isna().any().any() or frame[sample].duplicated().any():
        raise ValueError("QC sample IDs must be present and unique; isolate IDs must be present")
    fraction = pd.to_numeric(frame[mapping], errors="raise").astype(float)
    counts = pd.to_numeric(frame[assigned], errors="raise").astype(float)
    if mapping_total:
        total = pd.to_numeric(frame[mapping_total], errors="raise").astype(float)
        if not np.isfinite(total).all() or (total <= 0).any():
            raise ValueError("QC mapping denominators must be finite and positive")
        fraction = fraction / total
    if not np.isfinite(fraction).all() or not fraction.between(0, 1).all():
        raise ValueError("Mapping fraction must be in [0,1]; use --mapping-total-column for numerator counts")
    if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
        raise ValueError("Assigned-fragment counts must be nonnegative integers")
    data = pd.DataFrame({"sample": frame[sample].astype(str), "isolate": frame[isolate].astype(str),
                         "mapping_fraction": fraction, "assigned_fragments": counts})
    summary = data.groupby("isolate", sort=False)["assigned_fragments"].agg(["size", "mean", "std"])
    summary["CV_percent"] = np.where((summary["size"] >= 2) & (summary["mean"] > 0),
                                      100*summary["std"]/summary["mean"], np.nan)
    return data, summary


def qc_figure(frame, **columns):
    data, summary = qc_values(frame, **columns)
    palette = {name:plt.get_cmap("tab20")(i % 20) for i,name in enumerate(summary.index)}
    bar_colors = [palette[name] for name in data["isolate"]]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    x = np.arange(len(data))
    axes[0].bar(x, 100*data.mapping_fraction, color=bar_colors, width=.8)
    axes[0].set_ylabel("Mapping rate (%)"); axes[0].set_ylim(0, 100)
    axes[1].bar(x, data.assigned_fragments, color=bar_colors, width=.8)
    axes[1].set_ylabel("Assigned fragments")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for ax in axes[:2]:
        ax.set_xticks(x, data["sample"], rotation=90, fontsize=5)
    valid = summary["CV_percent"].notna()
    if valid.any():
        names = summary.index[valid]
        axes[2].bar(range(len(names)), summary.loc[names, "CV_percent"], color=[palette[name] for name in names])
        axes[2].set_xticks(range(len(names)), names, rotation=90, fontsize=6)
    else:
        axes[2].text(.5, .5, "CV needs at least two libraries\nwith a positive mean per isolate",
                     transform=axes[2].transAxes, ha="center", va="center", fontsize=7)
        axes[2].set_xticks([])
    axes[2].set_ylabel("Within-isolate CV (%)")
    for ax,title in zip(axes, ["a  Mapping", "b  Fragment assignment", "c  Replicate variation"]):
        ax.set_title(title, loc="left", fontweight="bold"); clean_axis(ax)
    return fig


def ratio_value(value):
    text = str(value).strip()
    if "/" in text:
        match = re.fullmatch(r"(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)", text)
        if not match or float(match[2]) <= 0:
            raise ValueError("Invalid gene ratio")
        result = float(match[1]) / float(match[2])
    else:
        result = float(text)
    if not np.isfinite(result) or not 0 <= result <= 1:
        raise ValueError("Gene ratios must be finite values in [0,1]")
    return result


def enrichment_figure(frame, *, top, term="term", ratio="gene_ratio", count="count", adjusted_p="padj", title="Supplied enrichment results"):
    needed = {term, ratio, count, adjusted_p}
    if frame.empty or not needed.issubset(frame) or top < 1:
        raise ValueError("Enrichment plotting requires its four columns and a positive explicit --top")
    if frame[term].isna().any():
        raise ValueError("Enrichment term names must be present")
    data = pd.DataFrame({"term": frame[term].astype(str), "ratio": frame[ratio].map(ratio_value),
                         "count": pd.to_numeric(frame[count], errors="raise"),
                         "padj": pd.to_numeric(frame[adjusted_p], errors="raise")})
    if not np.isfinite(data[["count", "padj"]]).all().all() or not data.padj.between(0, 1).all():
        raise ValueError("Supplied adjusted P values must be finite and in [0,1]")
    if (data["count"] <= 0).any() or not np.equal(data["count"], np.floor(data["count"])).all():
        raise ValueError("Gene counts must be positive integers")
    shown = data.sort_values("padj", kind="stable").head(top)
    positive = shown.loc[shown.padj > 0, "padj"]
    floor = max(float(positive.min())/10, np.finfo(float).tiny) if len(positive) else np.finfo(float).tiny
    logp = -np.log10(shown.padj.clip(lower=floor))
    max_count = float(shown["count"].max())
    fig, ax = plt.subplots(figsize=(7, max(2.8, .27*len(shown)+1.5)), layout="constrained")
    points = ax.scatter(shown.ratio, range(len(shown)), s=200*shown["count"]/max_count,
                         c=logp, cmap="viridis", norm=colors.Normalize(vmin=0, vmax=max(.1, float(logp.max()))),
                         edgecolor="#555555", linewidth=.35)
    ax.set_yticks(range(len(shown)), [textwrap.fill(value, 42) for value in shown.term])
    ax.invert_yaxis(); ax.set_xlabel("Gene ratio"); ax.set_title(title, loc="left")
    colorbar = fig.colorbar(points, ax=ax, fraction=.04, pad=.025)
    colorbar.set_label("-log10 supplied adjusted P")
    observed_counts = sorted(set(shown["count"].astype(int)))
    levels = [observed_counts[i] for i in sorted({0, len(observed_counts)//2, len(observed_counts)-1})]
    handles = [ax.scatter([], [], s=200*value/max_count, color="#999999", edgecolor="#555555", linewidth=.35) for value in levels]
    ax.legend(handles, list(map(str, levels)), title="Gene count", frameon=False,
               fontsize=6, title_fontsize=7, loc="upper left", bbox_to_anchor=(1.04, .45))
    ax.margins(x=.12, y=.12); clean_axis(ax)
    return fig


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.
    children: list = field(default_factory=list)


def read_tree(path):
    tokens = re.findall(r"'(?:[^']|'')*'|[(),:;]|[^\s(),:;]+", Path(path).read_text(encoding="utf-8"))
    position = 0
    def parse():
        nonlocal position
        if position >= len(tokens):
            raise ValueError("Truncated Newick tree")
        node = Node()
        if tokens[position] == "(":
            position += 1
            while True:
                node.children.append(parse())
                if position < len(tokens) and tokens[position] == ",":
                    position += 1
                elif position < len(tokens) and tokens[position] == ")":
                    position += 1; break
                else:
                    raise ValueError("Invalid Newick child list")
        if position < len(tokens) and tokens[position] not in {":", ",", ")", ";"}:
            node.name = tokens[position].strip("'").replace("''", "'")
            position += 1
        if position < len(tokens) and tokens[position] == ":":
            position += 1
            node.length = float(tokens[position]); position += 1
            if not np.isfinite(node.length) or node.length < 0:
                raise ValueError("Tree branch lengths must be finite and nonnegative")
        return node
    root = parse()
    if position != len(tokens)-1 or tokens[position] != ";":
        raise ValueError("A single Newick tree ending with a semicolon is required")
    return root


def tree_layout(root):
    tips, nodes = [], []
    def collect(node):
        nodes.append(node)
        if not node.children:
            tips.append(node)
        for child in node.children:
            collect(child)
    collect(root)
    labels = [node.name for node in tips]
    if not all(labels) or len(set(labels)) != len(labels):
        raise ValueError("Tree tip identifiers must be present and unique")
    branch_lengths = any(node.length > 0 for node in nodes)
    coords = {node: (0., float(i)) for i, node in enumerate(tips)}
    def place(node, x=0.):
        if node.children:
            children = [place(child, x+(child.length if branch_lengths else 1.)) for child in node.children]
            coords[node] = (x, (min(y for _,y in children)+max(y for _,y in children))/2)
        else:
            coords[node] = (x, coords[node][1])
        return coords[node]
    place(root)
    return labels, coords, branch_lengths


def mic_numbers(frame):
    def parse(value):
        if pd.isna(value):
            return np.nan
        match = re.fullmatch(r"\s*(?:<=|>=|[<>≤≥=])?\s*(\d+(?:\.\d+)?)\s*", str(value))
        if not match or float(match[1]) <= 0:
            raise ValueError(f"Invalid MIC display value {value!r}")
        return float(match[1])
    values = frame.map(parse).to_numpy(float)
    display = np.full(values.shape, np.nan)
    for j in range(values.shape[1]):
        good = np.isfinite(values[:, j])
        if not good.any():
            continue
        data = np.log2(values[good, j]); span = data.max()-data.min()
        display[good, j] = (data-data.min())/span if span > 0 else .5
    return display


def integrated_figure(tree, metadata, mics, states=None, qmatrices=(), *, fields=("source","location","date")):
    labels, coords, branch_lengths = tree_layout(tree)
    if not set(labels).issubset(metadata.index):
        raise ValueError("Metadata must include every tree tip")
    fields = [name for name in fields if name in metadata]
    if not fields:
        raise ValueError("No requested metadata columns are present")
    if not set(mics.index).issubset(labels):
        raise ValueError("MIC table contains samples outside the supplied tree")
    if states is not None and not set(states.index).issubset(labels):
        raise ValueError("State table contains samples outside the supplied tree")
    widths = [3.1, max(1.5, .85*len(fields)), max(1.4, .22*len(mics.columns))]
    if states is not None:
        widths.append(max(.8, .25*len(states.columns)))
    widths.extend([.65]*len(qmatrices))
    fig = plt.figure(figsize=(sum(widths)+1.3, max(3, len(labels)*.17+1.5)), layout="constrained")
    grid = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=.04)
    axes = [fig.add_subplot(grid[0,i]) for i in range(len(widths))]
    ax = axes[0]
    span = max(x for x,y in coords.values()) or 1.
    for node,(x,y) in coords.items():
        if node.children:
            child_y = [coords[c][1] for c in node.children]
            ax.plot([x,x], [min(child_y),max(child_y)], color="black", lw=.6)
            for child in node.children:
                cx,cy = coords[child]
                ax.plot([x,cx],[cy,cy],color="black",lw=.6)
                if child.children and child.name:
                    ax.text(cx,cy-.12,child.name,fontsize=4.5,ha="right",va="bottom")
    for name,node in zip(labels,[node for node in coords if not node.children]):
        x,y = coords[node]
        ax.text(x+span*.025,y,name,fontsize=6,va="center")
    ax.set_xlim(-span*.02,span*1.42); ax.set_yticks([])
    ax.set_xlabel("Substitutions per site" if branch_lengths else "Topology depth")
    ax.set_title("a  Phylogeny",loc="left",fontweight="bold")
    clean_axis(ax); ax.spines['left'].set_visible(False)
    meta = metadata.loc[labels]
    ax = axes[1]
    for j,name in enumerate(fields):
        for i,value in enumerate(meta[name]):
            ax.text(j,i,"—" if pd.isna(value) else str(value),ha="center",va="center",fontsize=5.5)
    ax.set_xticks(range(len(fields)),fields,rotation=45,ha="right")
    ax.set_xlim(-.55,len(fields)-.45); ax.set_yticks([])
    ax.set_title("b  Metadata",loc="left",fontweight="bold")
    for spine in ax.spines.values():spine.set_visible(False)
    ax = axes[2]
    palette=plt.get_cmap('Blues').copy();palette.set_bad(MISSING)
    ax.imshow(mic_numbers(mics.reindex(labels)),vmin=0,vmax=1,cmap=palette,aspect="auto",interpolation="nearest")
    ax.set_xticks(range(len(mics.columns)),mics.columns,rotation=90);ax.set_yticks([])
    ax.set_title("c  MICs",loc="left",fontweight="bold")
    cursor=3
    if states is not None:
        ax=axes[cursor];cursor+=1
        shown=states.reindex(labels)
        for i in range(len(labels)):
            for j in range(len(shown.columns)):
                value=shown.iloc[i,j]
                if pd.isna(value):
                    ax.plot(j,i,marker='x',color='#777777',ms=3,mew=.6)
                elif str(value) in {'0','0.0','1','1.0'}:
                    ax.plot(j,i,marker='o',mfc='black' if float(value)==1 else 'white',mec='black',ms=3,mew=.6)
                else:raise ValueError('Amino-acid state values must be 0 (reference), 1 (variant), or missing')
        ax.set_xlim(-.6,len(shown.columns)-.4);ax.set_yticks([])
        ax.set_xticks(range(len(shown.columns)),shown.columns,rotation=90)
        ax.set_title("d  AA states",loc="left",fontweight="bold")
    for number,q in enumerate(qmatrices):
        if not set(labels).issubset(q.index):raise ValueError('Q matrix must include every tree tip')
        values=q.loc[labels].to_numpy(float)
        if not np.isfinite(values).all() or (values<0).any() or (values>1).any() or not np.allclose(values.sum(axis=1),1,atol=.002):
            raise ValueError('Q matrix must contain finite proportions summing to one per sample')
        values=values/values.sum(axis=1,keepdims=True)
        ax=axes[cursor];cursor+=1;left=np.zeros(len(labels))
        for k in range(values.shape[1]):
            ax.barh(np.arange(len(labels)),values[:,k],left=left,height=.95,
                    color=plt.get_cmap('tab10')(k%10),linewidth=0)
            left+=values[:,k]
        ax.set_xlim(0,1);ax.set_xticks([0,1]);ax.set_yticks([])
        ax.set_title(('e  ' if number==0 else '')+f'K={values.shape[1]}',fontsize=7)
    for ax in axes:
        ax.set_ylim(len(labels)-.5,-.5);ax.tick_params(axis='y',length=0)
    return fig


def main():
    configure()
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    def output(p):
        p.add_argument('--output',type=Path,required=True);p.add_argument('--dpi',type=int,default=600)
    p=sub.add_parser('ani');p.add_argument('--matrix',type=Path,required=True);p.add_argument('--subset',type=Path)
    p.add_argument('--vmin',type=float,default=86.);p.add_argument('--vmax',type=float,default=100.);output(p)
    p=sub.add_parser('rna');p.add_argument('--coordinates',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--variance',type=Path);p.add_argument('--differential-expression',type=Path,required=True)
    p.add_argument('--gene-sets',type=Path,required=True);p.add_argument('--effect-column',default='log2FoldChange')
    p.add_argument('--adjusted-p-column',default='padj');p.add_argument('--flip-pc1',action='store_true');output(p)
    p=sub.add_parser('rms');p.add_argument('--matrix',type=Path,required=True);output(p)
    p=sub.add_parser('heatmap');p.add_argument('--matrix',type=Path,required=True);p.add_argument('--genes',type=Path)
    p.add_argument('--log2',action='store_true');p.add_argument('--row-zscore',action='store_true');output(p)
    p=sub.add_parser('integrated');p.add_argument('--tree',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--mic',type=Path,required=True);p.add_argument('--states',type=Path);p.add_argument('--qmatrix',type=Path,action='append',default=[])
    p.add_argument('--fields',nargs='+',default=['source','location','date']);output(p)
    p=sub.add_parser('variants');p.add_argument('--table',type=Path,required=True)
    p.add_argument('--contig-column',default='CHROM');p.add_argument('--position-column',default='POS')
    p.add_argument('--impact-column',default='impact');p.add_argument('--gene-column',default='gene_id')
    p.add_argument('--label-genes',action='store_true');output(p)
    p=sub.add_parser('qc');p.add_argument('--table',type=Path,required=True)
    p.add_argument('--sample-column',default='unit');p.add_argument('--isolate-column',default='isolate')
    p.add_argument('--mapping-column',default='mapped_fraction');p.add_argument('--assigned-column',default='assigned_fragments')
    p.add_argument('--mapping-total-column');output(p)
    p=sub.add_parser('enrichment');p.add_argument('--table',type=Path,required=True);p.add_argument('--top',type=int,required=True)
    p.add_argument('--term-column',default='term');p.add_argument('--ratio-column',default='gene_ratio')
    p.add_argument('--count-column',default='count');p.add_argument('--adjusted-p-column',default='padj')
    p.add_argument('--title',default='Supplied enrichment results');output(p)
    args=parser.parse_args()
    if args.dpi<72:raise ValueError('DPI must be at least 72')
    if args.command=='ani':
        fig=ani_figure(table(args.matrix,numeric=True),ids_file(args.subset) if args.subset else None,vmin=args.vmin,vmax=args.vmax)
    elif args.command=='rna':
        fig=rna_figure(table(args.coordinates,numeric=True),table(args.metadata),
            pd.read_csv(args.differential_expression),pd.read_csv(args.gene_sets),
            pd.read_csv(args.variance) if args.variance else None,flip_pc1=args.flip_pc1,
            effect=args.effect_column,adjusted_p=args.adjusted_p_column)
    elif args.command=='rms':
        frame=square_matrix(table(args.matrix,numeric=True),symmetry=True)
        if (frame.to_numpy()<0).any():raise ValueError('RMS distances cannot be negative')
        fig,ax=plt.subplots(figsize=(5,4.5),layout='constrained')
        im=matrix_panel(ax,frame,'Expression RMS distance',vmin=0,vmax=float(frame.max().max()) or 1.,annotate=True)
        fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
    elif args.command=='heatmap':
        fig=heatmap_figure(table(args.matrix,numeric=True),genes=ids_file(args.genes) if args.genes else None,
            log2=args.log2,row_zscore=args.row_zscore)
    elif args.command=='variants':
        fig=variants_figure(pd.read_csv(args.table,sep=None,engine='python'),contig=args.contig_column,
            position=args.position_column,impact=args.impact_column,gene=args.gene_column,label_genes=args.label_genes)
    elif args.command=='qc':
        fig=qc_figure(pd.read_csv(args.table,sep=None,engine='python'),sample=args.sample_column,
            isolate=args.isolate_column,mapping=args.mapping_column,assigned=args.assigned_column,
            mapping_total=args.mapping_total_column)
    elif args.command=='enrichment':
        fig=enrichment_figure(pd.read_csv(args.table,sep=None,engine='python'),top=args.top,term=args.term_column,
            ratio=args.ratio_column,count=args.count_column,adjusted_p=args.adjusted_p_column,title=args.title)
    else:
        fig=integrated_figure(read_tree(args.tree),table(args.metadata),table(args.mic),
            table(args.states) if args.states else None,[table(p,numeric=True) for p in args.qmatrix],fields=args.fields)
    export(fig,args.output,args.dpi)
    print(args.output)


if __name__=='__main__':
    main()
