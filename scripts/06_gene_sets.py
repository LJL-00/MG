#!/usr/bin/env python3
"""Analyse predefined gene sets by moderated t statistics and exact label permutations.

Expression CSV uses Gene ID, optional Gene Name/Gene Description, and replicate
TPM columns named <isolate>_<replicate number>. Groups CSV uses isolate and group
(environment or clinical). Expression filtering requires mean TPM >= 1 in at
least two isolates. Replicate log2(TPM+1) values are averaged within isolate.
All distinct group assignments are enumerated; four versus four gives 70.
"""
from __future__ import annotations
from itertools import combinations
from math import log
from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd

STRAINS = []
DISPLAY = {}
ENVIRONMENT = set()
CLINICAL = set()
EXPR_PATH = None
MEMBERSHIP_PATH = None

MAIN_ORDER = [
    "Ergosterol_biosynthesis",
    "ABC_efflux",
    "MFS_xenobiotic_efflux",
    "Chitin_synthesis",
    "Mannan_glycosylation",
    "ROS_detoxification",
    "HOG_MAPK_osmotic",
    "Heat_shock_chaperone",
    "Phosphate_nutrient_assimilation",
    "ALS_adhesion",
    "SAP_aspartic_protease",
    "Hyphal_polarity_calcium",
]
SUPPLEMENTAL_ORDER = [
    "Broad_cell_wall_integrity",
    "V_ATPase_vacuolar_acidification",
]
ALL_MODULE_ORDER = [
    "Ergosterol_biosynthesis",
    "ABC_efflux",
    "MFS_xenobiotic_efflux",
    "Chitin_synthesis",
    "Mannan_glycosylation",
    "Broad_cell_wall_integrity",
    "V_ATPase_vacuolar_acidification",
    "ROS_detoxification",
    "HOG_MAPK_osmotic",
    "Heat_shock_chaperone",
    "Phosphate_nutrient_assimilation",
    "ALS_adhesion",
    "SAP_aspartic_protease",
    "Hyphal_polarity_calcium",
]


def digamma(x: float) -> float:
    result = 0.0
    while x < 8.0:
        result -= 1.0 / x
        x += 1.0
    inv = 1.0 / x
    inv2 = inv * inv
    result += (
        log(x)
        - 0.5 * inv
        - inv2 / 12.0
        + inv2 * inv2 / 120.0
        - inv2 * inv2 * inv2 / 252.0
    )
    return result


def trigamma(x: float) -> float:
    result = 0.0
    while x < 8.0:
        result += 1.0 / (x * x)
        x += 1.0
    inv = 1.0 / x
    inv2 = inv * inv
    result += (
        inv
        + 0.5 * inv2
        + inv2 * inv / 6.0
        - inv2 * inv2 * inv / 30.0
        + inv2 * inv2 * inv2 * inv / 42.0
    )
    return result


def inverse_trigamma(value: float) -> float:
    if value <= 0:
        return 5e5
    lo, hi = 1e-4, 5e5
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if trigamma(mid) > value:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def bh_fdr(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    p = values[finite]
    if len(p) == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    out[finite] = restored
    return out


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2.0 + 1.0
        ranks[order[i:j]] = rank
        i = j
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = average_ranks(np.asarray(x, dtype=float))
    ry = average_ranks(np.asarray(y, dtype=float))
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def fit_variance_prior(s2: np.ndarray, df: int) -> tuple[float, float]:
    finite = s2[np.isfinite(s2) & (s2 > 1e-12)]
    if len(finite) < 2:
        raise ValueError("At least two genes with positive residual variance are needed")
    log_s2 = np.log(finite)
    residual = max(float(np.var(log_s2, ddof=1)) - trigamma(df / 2.0), 1e-8)
    d0 = min(1e6, 2.0 * inverse_trigamma(residual))
    log_s0 = (
        float(np.mean(log_s2))
        - digamma(df / 2.0)
        + log(df / 2.0)
        + digamma(d0 / 2.0)
        - log(d0 / 2.0)
    )
    return float(np.exp(log_s0)), float(d0)


def moderated_statistics(
    matrix: np.ndarray,
    clinical_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    clinical = matrix[:, clinical_mask]
    environment = matrix[:, ~clinical_mask]
    n_clinical = clinical.shape[1]
    n_environment = environment.shape[1]
    effect = clinical.mean(axis=1) - environment.mean(axis=1)
    df = n_clinical + n_environment - 2
    s2 = (
        (n_clinical - 1) * clinical.var(axis=1, ddof=1)
        + (n_environment - 1) * environment.var(axis=1, ddof=1)
    ) / df
    s0, d0 = fit_variance_prior(s2, df)
    posterior = (d0 * s0 + df * s2) / (d0 + df)
    se = np.sqrt(posterior * (1.0 / n_clinical + 1.0 / n_environment))
    tstat = np.divide(effect, se, out=np.zeros_like(effect), where=se > 0)
    return effect, tstat, s0, d0


def all_label_masks(n: int, n_clinical: int) -> list[np.ndarray]:
    masks = []
    for indices in combinations(range(n), n_clinical):
        mask = np.zeros(n, dtype=bool)
        mask[list(indices)] = True
        masks.append(mask)
    return masks


def enrichment_score(
    scores: np.ndarray,
    member_mask: np.ndarray,
) -> float:
    order = np.argsort(-scores)
    ordered_hits = member_mask[order]
    n_hits = int(ordered_hits.sum())
    if n_hits < 1 or n_hits == len(scores):
        return np.nan
    weights = np.abs(scores[order])
    hit_weights = weights * ordered_hits
    if hit_weights.sum() == 0:
        hit_weights = ordered_hits.astype(float)
    running = np.cumsum(
        hit_weights / hit_weights.sum()
        - (~ordered_hits).astype(float) / (len(scores) - n_hits)
    )
    max_value = float(running.max())
    min_value = float(running.min())
    return max_value if abs(max_value) >= abs(min_value) else min_value


def normalize_gene_set_name(value: str) -> str:
    replacements = {
        "Ergosterol_steroid_biosynthesis": "Ergosterol_biosynthesis",
        "Efflux_ABC_transporters": "ABC_efflux",
        "Efflux_MFS_xenobiotic": "MFS_xenobiotic_efflux",
        "Cell_wall_chitin": "Chitin_synthesis",
        "Cell_wall_mannan_glycosylation": "Mannan_glycosylation",
        "Stress_HOG_MAPK_osmotic": "HOG_MAPK_osmotic",
        "Stress_ROS_detoxification": "ROS_detoxification",
        "Stress_heat_shock_chaperone": "Heat_shock_chaperone",
        "Nutrient_phosphate_assimilation": "Phosphate_nutrient_assimilation",
        "Virulence_SAP_aspartic_protease": "SAP_aspartic_protease",
        "Virulence_ALS_adhesion": "ALS_adhesion",
        "Virulence_hyphal_polarity_calcium": "Hyphal_polarity_calcium",
    }
    return replacements.get(str(value), str(value))


def prepare_expression() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    expr = pd.read_csv(EXPR_PATH, encoding="utf-8-sig")
    if expr["Gene ID"].isna().any() or expr["Gene ID"].astype(str).duplicated().any():
        raise ValueError("Gene IDs must be present and unique")
    for column in ("Gene Name", "Gene Description"):
        if column not in expr:
            expr[column] = ""
    replicate_log, strain_log, strain_tpm = [], [], []
    for strain in STRAINS:
        pattern = re.compile(re.escape(strain) + r"_(\d+)$")
        matches = [(int(match[1]), column) for column in expr.columns
                   if (match := pattern.fullmatch(column))]
        columns = [column for _, column in sorted(matches)]
        if len(columns) < 2:
            raise ValueError(f"At least two replicate TPM columns are needed for {strain}")
        tpm = expr[columns].to_numpy(float)
        if not np.isfinite(tpm).all() or (tpm < 0).any():
            raise ValueError("TPM values must be finite and nonnegative")
        replicate_log.append(np.log2(tpm + 1.0))
        strain_log.append(np.log2(tpm + 1.0).mean(axis=1))
        strain_tpm.append(tpm.mean(axis=1))
    replicate_matrix = np.concatenate(replicate_log, axis=1)
    strain_matrix = np.column_stack(strain_log)
    strain_tpm_matrix = np.column_stack(strain_tpm)
    keep = (strain_tpm_matrix >= 1.0).sum(axis=1) >= 2
    if keep.sum() < 2:
        raise ValueError("Fewer than two genes pass the mean-TPM filter")
    meta = expr.loc[keep, ["Gene ID", "Gene Name", "Gene Description"]].reset_index(drop=True)
    return meta, strain_matrix[keep], replicate_matrix[keep]


def prepare_gene_sets(
    gene_ids: np.ndarray,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    membership = pd.read_csv(MEMBERSHIP_PATH, encoding="utf-8-sig")
    membership["gene_set"] = membership["gene_set"].map(normalize_gene_set_name)
    gene_index = {str(gene): i for i, gene in enumerate(gene_ids)}
    masks: dict[str, np.ndarray] = {}
    rows = []
    for gene_set, sub in membership.groupby("gene_set"):
        if gene_set not in set(ALL_MODULE_ORDER):
            continue
        mask = np.zeros(len(gene_ids), dtype=bool)
        for gene in sub["Gene ID"].astype(str):
            if gene in gene_index:
                mask[gene_index[gene]] = True
        masks[gene_set] = mask
        rows.append(
            {
                "gene_set": gene_set,
                "tier": sub["tier"].iloc[0] if "tier" in sub else "",
                "n_genes_in_filtered_matrix": int(mask.sum()),
            }
        )
    missing = set(MAIN_ORDER) - set(masks)
    if missing:
        raise ValueError(f"Missing primary gene sets: {sorted(missing)}")
    if any(not masks[name].any() or masks[name].all() for name in MAIN_ORDER):
        raise ValueError("Each primary set must contain some, but not all, filtered genes")
    return masks, pd.DataFrame(rows)


def module_scores(
    matrix: np.ndarray,
    masks: dict[str, np.ndarray],
) -> pd.DataFrame:
    mean = matrix.mean(axis=1, keepdims=True)
    sd = matrix.std(axis=1, ddof=0, keepdims=True)
    z = np.divide(matrix - mean, sd, out=np.zeros_like(matrix), where=sd > 0)
    rows = []
    for gene_set in ALL_MODULE_ORDER:
        if gene_set not in masks:
            continue
        mask = masks[gene_set]
        scores = z[mask].mean(axis=0)
        for index, strain in enumerate(STRAINS):
            rows.append(
                {
                    "gene_set": gene_set,
                    "strain": DISPLAY[strain],
                    "group": (
                        "environment"
                        if strain in ENVIRONMENT
                        else "clinical"
                    ),
                    "score": float(scores[index]),
                }
            )
    return pd.DataFrame(rows)


def analyse_scope(
    name: str,
    matrix: np.ndarray,
    strains: list[str],
    masks: dict[str, np.ndarray],
    score_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    clinical_mask = np.array([strain in CLINICAL for strain in strains])
    permutations = all_label_masks(len(strains), int(clinical_mask.sum()))
    observed_index = next(
        i for i, mask in enumerate(permutations) if np.array_equal(mask, clinical_mask)
    )
    permutation_effects = []
    permutation_t = []
    for mask in permutations:
        effect, tstat, _, _ = moderated_statistics(matrix, mask)
        permutation_effects.append(effect)
        permutation_t.append(tstat)
    permutation_effects = np.asarray(permutation_effects)
    permutation_t = np.asarray(permutation_t)
    observed_effect = permutation_effects[observed_index]
    observed_t = permutation_t[observed_index]

    score_pivot = score_table.pivot(
        index="gene_set", columns="strain", values="score"
    )
    display_strains = [DISPLAY[strain] for strain in strains]
    score_matrix = score_pivot.reindex(index=ALL_MODULE_ORDER, columns=display_strains).to_numpy(float)
    score_permutation_effect = np.asarray(
        [
            score_matrix[:, mask].mean(axis=1)
            - score_matrix[:, ~mask].mean(axis=1)
            for mask in permutations
        ]
    )

    rows = []
    gsea_rows = []
    for gene_set in MAIN_ORDER:
        member_mask = masks[gene_set]
        gene_effect_perm = permutation_effects[:, member_mask]
        positive_fraction_perm = (gene_effect_perm > 0).mean(axis=1)
        es_perm = np.asarray(
            [
                enrichment_score(permutation_t[i], member_mask)
                for i in range(len(permutations))
            ]
        )
        observed_es = es_perm[observed_index]
        same_sign = es_perm[es_perm >= 0] if observed_es >= 0 else es_perm[es_perm < 0]
        nes = (
            observed_es / np.mean(np.abs(same_sign))
            if len(same_sign) and np.mean(np.abs(same_sign)) > 0
            else np.nan
        )
        rows.append(
            {
                "analysis_scope": name,
                "gene_set": gene_set,
                "n_genes": int(member_mask.sum()),
                "module_score_effect": float(
                    score_permutation_effect[
                        observed_index, ALL_MODULE_ORDER.index(gene_set)
                    ]
                ),
                "module_score_permutation_p": float(
                    (
                        np.abs(
                            score_permutation_effect[:, ALL_MODULE_ORDER.index(gene_set)]
                        )
                        >= abs(
                            score_permutation_effect[
                                observed_index, ALL_MODULE_ORDER.index(gene_set)
                            ]
                        )
                        - 1e-12
                    ).mean()
                ),
                "positive_fraction": float(
                    positive_fraction_perm[observed_index]
                ),
                "directionality_permutation_p": float(
                    (
                        np.abs(positive_fraction_perm - 0.5)
                        >= abs(positive_fraction_perm[observed_index] - 0.5)
                        - 1e-12
                    ).mean()
                ),
                "moderated_t_mean": float(observed_t[member_mask].mean()),
                "GSEA_ES": float(observed_es),
                "GSEA_NES": float(nes),
                "GSEA_phenotype_permutation_p": float(
                    (np.abs(es_perm) >= abs(observed_es) - 1e-12).mean()
                ),
                "n_label_permutations": len(permutations),
            }
        )
        gsea_rows.extend(
            {
                "analysis_scope": name,
                "gene_set": gene_set,
                "permutation_index": index,
                "is_observed": index == observed_index,
                "ES": float(value),
            }
            for index, value in enumerate(es_perm)
        )

    for gene_set in SUPPLEMENTAL_ORDER:
        if gene_set not in masks:
            continue
        index = ALL_MODULE_ORDER.index(gene_set)
        observed = score_permutation_effect[observed_index, index]
        rows.append(
            {
                "analysis_scope": name,
                "gene_set": gene_set,
                "n_genes": int(masks[gene_set].sum()),
                "module_score_effect": float(observed),
                "module_score_permutation_p": float(
                    (
                        np.abs(score_permutation_effect[:, index])
                        >= abs(observed) - 1e-12
                    ).mean()
                ),
                "positive_fraction": np.nan,
                "directionality_permutation_p": np.nan,
                "moderated_t_mean": np.nan,
                "GSEA_ES": np.nan,
                "GSEA_NES": np.nan,
                "GSEA_phenotype_permutation_p": np.nan,
                "n_label_permutations": len(permutations),
            }
        )

    results = pd.DataFrame(rows)
    main_mask = results["gene_set"].isin(MAIN_ORDER)
    results.loc[main_mask, "module_score_fdr"] = bh_fdr(
        results.loc[main_mask, "module_score_permutation_p"].to_numpy(float)
    )
    results.loc[main_mask, "directionality_fdr"] = bh_fdr(
        results.loc[main_mask, "directionality_permutation_p"].to_numpy(float)
    )
    results.loc[main_mask, "GSEA_fdr"] = bh_fdr(
        results.loc[main_mask, "GSEA_phenotype_permutation_p"].to_numpy(float)
    )
    return results, pd.DataFrame(gsea_rows), observed_t


def main() -> None:
    global STRAINS, DISPLAY, ENVIRONMENT, CLINICAL, EXPR_PATH, MEMBERSHIP_PATH
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    groups = pd.read_csv(args.groups, dtype=str)
    if groups[["isolate", "group"]].isna().any().any() or groups.isolate.duplicated().any():
        raise ValueError("Each isolate must have one unique group assignment")
    if not set(groups.group) == {"environment", "clinical"}:
        raise ValueError("Expected environment and clinical groups")
    if groups.group.value_counts().min() < 2:
        raise ValueError("At least two isolates per group are required")
    STRAINS = groups.isolate.tolist()
    DISPLAY = {name: name for name in STRAINS}
    ENVIRONMENT = set(groups.loc[groups.group == "environment", "isolate"])
    CLINICAL = set(groups.loc[groups.group == "clinical", "isolate"])
    EXPR_PATH, MEMBERSHIP_PATH = args.expression, args.membership
    meta, strain_matrix, _ = prepare_expression()
    masks, membership_summary = prepare_gene_sets(meta["Gene ID"].astype(str).to_numpy())
    scores = module_scores(strain_matrix, masks)
    results, nulls, tstat = analyse_scope("all_isolates", strain_matrix, STRAINS, masks, scores)
    results = results[results["gene_set"].isin(MAIN_ORDER)].copy()
    args.output.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output / "gene_set_summaries.csv", index=False)
    rank = meta.copy()
    rank["moderated_t"] = tstat
    rank.sort_values("moderated_t", ascending=False).to_csv(args.output / "moderated_t_ranked_genes.csv", index=False)
    membership_summary.to_csv(args.output / "gene_set_membership_summary.csv", index=False)
    nulls.to_csv(args.output / "exact_label_permutation_nulls.csv", index=False)
    print(f"Analysed {len(MAIN_ORDER)} gene sets in {len(STRAINS)} isolates")


if __name__ == "__main__":
    main()
