"""Auditable local sample QC, statistical tests and gene-set over-representation.

PCA uses centered, unscaled log2 protein abundances and NumPy SVD. Missing-value
filling is optional and applies to PCA only, never to a differential test.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import re
import warnings as pywarnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy import stats
from plotly.offline import get_plotlyjs

from date_display import display_provenance
from figure_style import COLORS, CATEGORICAL, DIVERGING, publication_layout, math_markup, math_assets
from proteomics import _read, _number, _clean, _text, _detect_mapping, _qc_columns, _truthy_qc, _csv_safe, bh_adjust


def run_test(control, treatment, method="welch"):
    """Two-sided test; t statistics and mean differences are treatment−control.

    Values are on the caller's chosen analysis scale. Nonfinite independent
    observations are omitted. Paired observations are matched by their original
    positions and omitted jointly. Unestimable tests return p_value=None with a
    note, while an invalid method or unequal paired lengths raises ValueError.
    """
    aliases = {"ttest": "student", "student_t": "student", "paired_t": "paired",
               "mann_whitney": "mannwhitney", "mann-whitney": "mannwhitney"}
    method = aliases.get(method, method)
    if method not in {"welch", "student", "paired", "mannwhitney"}:
        raise ValueError("Test must be welch, student, paired or mannwhitney.")
    def numeric(values):
        try:
            original = list(values)
        except TypeError as exc:
            raise ValueError("Test inputs must be one-dimensional sample arrays.") from exc
        if any(isinstance(value, (list, tuple, np.ndarray)) for value in original):
            raise ValueError("Test inputs must be one-dimensional sample arrays.")
        return np.asarray([np.nan if _number(v) is None else _number(v) for v in original], dtype=float)
    a, b = numeric(control), numeric(treatment)
    original_counts = (len(a), len(b))
    if method == "paired":
        if len(a) != len(b):
            raise ValueError("Paired t tests require equal group sizes, paired in column-list order.")
        keep = np.isfinite(a) & np.isfinite(b)
        a, b = a[keep], b[keep]
    else:
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    result = {"method": method, "statistic": None, "p_value": None,
              "n_control": len(a), "n_treatment": len(b),
              "n_pairs": len(a) if method == "paired" else None,
              "missing_control": original_counts[0] - len(a), "missing_treatment": original_counts[1] - len(b),
              "mean_control": float(np.mean(a)) if len(a) else None,
              "mean_treatment": float(np.mean(b)) if len(b) else None, "mean_difference": None, "note": ""}
    if len(a) and len(b):
        result["mean_difference"] = float(np.mean(b) - np.mean(a))
    if min(len(a), len(b)) < 2:
        result["note"] = "Fewer than 2 complete pairs; p not calculated" if method == "paired" else "Fewer than 2 valid replicates per group; p not calculated"
        return result
    if method == "paired" and np.var(b - a, ddof=1) == 0:
        result["note"] = "Zero variance in paired differences; t test not estimable"
        return result
    if method in {"welch", "student"} and np.var(a, ddof=1) == 0 and np.var(b, ddof=1) == 0:
        result["note"] = "Zero variance in both groups; t test not estimable"
        return result
    with pywarnings.catch_warnings():
        pywarnings.simplefilter("ignore", RuntimeWarning)
        if method == "paired":
            test = stats.ttest_rel(b, a, alternative="two-sided")
        elif method == "mannwhitney":
            test = stats.mannwhitneyu(b, a, alternative="two-sided", method="auto")
        else:
            test = stats.ttest_ind(b, a, equal_var=method == "student", alternative="two-sided")
    result.update(statistic=_clean(test.statistic), p_value=_clean(test.pvalue))
    if result["p_value"] is None:
        result["note"] = "The test returned no finite p value"
    if method == "mannwhitney":
        result["note"] = "Two-sided rank test, not a test of mean differences. SciPy auto uses an exact method for small samples without ties; otherwise asymptotic."
    return result


def _matrix(frame, config):
    mapping = {**_detect_mapping(frame.columns), **(config.get("mapping") or {})}
    control, treatment = config.get("control") or [], config.get("treatment") or []
    if not isinstance(control, list) or not isinstance(treatment, list):
        raise ValueError("Replicate columns must be lists.")
    sample_groups = config.get("sample_groups")
    columns = (mapping.get("abundance") or []) if sample_groups is not None else (control + treatment if control or treatment else mapping.get("abundance") or [])
    if not isinstance(columns, list) or not columns or len(columns) != len(set(columns)):
        raise ValueError("Select at least one distinct abundance column. Groups cannot share columns.")
    if any(col not in frame.columns for col in columns):
        raise ValueError("Abundance column not found.")
    if sample_groups is not None:
        if not isinstance(sample_groups, dict) or set(sample_groups) != set(columns):
            raise ValueError("sample_groups must assign every selected abundance column exactly once.")
        if any(not isinstance(group, str) or not group.strip() for group in sample_groups.values()):
            raise ValueError("Each sample group must have a nonempty label.")
    scale, normalization = config.get("abundance_scale", "raw"), config.get("normalization", "none")
    if scale not in {"raw", "log2", "log10", "ln"}:
        raise ValueError("Abundance scale must be raw, log2, log10 or ln.")
    if normalization not in {"none", "median"}:
        raise ValueError("Normalization must be none or median.")
    values = [[_number(v) for v in frame[col]] for col in columns]
    matrix = np.array(values, dtype=float).T
    if scale == "raw":
        matrix[matrix <= 0] = np.nan
        with np.errstate(invalid="ignore", divide="ignore"):
            matrix = np.log2(matrix)
    elif scale == "log10":
        matrix *= math.log2(10)
    elif scale == "ln":
        matrix /= math.log(2)
    flags = _qc_columns(frame.columns)
    excluded = np.array([bool(config.get("exclude_qc") and any(_truthy_qc(row.get(c)) for c in flags))
                         for row in frame.to_dict(orient="records")], dtype=bool)
    offsets = {}
    if normalization == "median":
        with pywarnings.catch_warnings():
            pywarnings.simplefilter("ignore", RuntimeWarning)
            medians = np.nanmedian(matrix[~excluded], axis=0)
        if not np.isfinite(medians).all():
            raise ValueError("A sample has no valid abundance for median normalization.")
        adjustment = medians - np.median(medians)
        matrix -= adjustment
        offsets = dict(zip(columns, map(float, adjustment)))
    groups = ([sample_groups[c].strip() for c in columns] if sample_groups is not None else
              ["control" if c in control else "treatment" if c in treatment else "ungrouped" for c in columns])
    return matrix, columns, groups, excluded, offsets, mapping


def _pca(matrix, row_ids, policy):
    if policy not in {"complete", "mean"}:
        raise ValueError("PCA missing-data policy must be complete or mean.")
    eligible = np.isfinite(matrix).all(axis=1) if policy == "complete" else np.isfinite(matrix).any(axis=1)
    selected = matrix[eligible].copy()
    feature_ids = np.asarray(row_ids)[eligible].tolist()
    missing_count = int((~np.isfinite(selected)).sum())
    result = {"available": False, "policy": policy, "feature_count": len(selected),
              "imputed_values": missing_count, "feature_row_ids": feature_ids,
              "scores": [], "loadings": [], "explained_variance": [], "explained_variance_ratio": [], "note": ""}
    if matrix.shape[1] < 2 or len(selected) < 2:
        result["note"] = "PCA needs at least 2 samples and 2 eligible protein rows."
        return result
    if policy == "mean":
        means = np.nanmean(selected, axis=1)
        positions = np.where(~np.isfinite(selected))
        selected[positions] = means[positions[0]]
    centered = selected.T - np.mean(selected.T, axis=0)
    if not np.any(centered):
        result["note"] = "PCA is unavailable: selected proteins have no between-sample variance."
        return result
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    components = min(2, len(singular))
    # Choose a stable orientation: the largest absolute loading is positive.
    for index in range(components):
        pivot = int(np.argmax(np.abs(vt[index])))
        if vt[index, pivot] < 0:
            u[:, index] *= -1
            vt[index] *= -1
    variance = singular ** 2 / (matrix.shape[1] - 1)
    scores = u[:, :components] * singular[:components]
    rank = int(np.sum(singular > singular[0] * max(centered.shape) * np.finfo(float).eps))
    result.update(available=True, scores=scores.tolist(), loadings=vt[:components].T.tolist(),
                  explained_variance=variance[:components].tolist(),
                  explained_variance_ratio=(variance[:components] / variance.sum()).tolist(), rank=rank)
    if rank < 2:
        result["note"] = "Only one nonzero component; PC2 carries no independent variation."
    return result


def _raw_cv(log2_values):
    """Sample CV on reconstructed linear abundance, stabilized against overflow."""
    values = np.asarray(log2_values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return None
    scaled = np.exp2(values - np.max(values))
    return float(np.std(scaled, ddof=1) / np.mean(scaled))


def _write_table(path, rows, columns):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: _csv_safe(row.get(key)) for key in columns} for row in rows)


def _layout(title, **kwargs):
    return publication_layout(title, **kwargs)


def analyze_quality(path, config, output_dir):
    """Read selected abundance columns and write a self-contained offline QC report."""
    config = config or {}
    frame, _, chosen, header = _read(path, config.get("sheet"), config.get("header"))
    matrix, columns, groups, excluded, offsets, mapping = _matrix(frame, config)
    notices = []
    if "ungrouped" in groups:
        notices.append("Ungrouped CV measures variation across selected samples, not technical replicate variation.")
    if excluded.any():
        notices.append(f"{int(excluded.sum())} QC-flagged rows excluded from QC and PCA.")
    selected = matrix[~excluded]
    total = len(selected)
    row_ids = np.arange(header + 2, header + 2 + len(frame))[~excluded]
    metrics = []
    for index, col in enumerate(columns):
        values = selected[:, index]
        finite = values[np.isfinite(values)]
        metrics.append({"sample": col, "group": groups[index], "eligible_protein_rows": total,
                        "observed": len(finite), "missing": total - len(finite),
                        "missing_fraction": (total - len(finite)) / total if total else None,
                        "median_log2": float(np.median(finite)) if len(finite) else None,
                        "mean_log2": float(np.mean(finite)) if len(finite) else None,
                        "sd_log2": float(np.std(finite, ddof=1)) if len(finite) > 1 else None,
                        "normalization_offset_log2": offsets.get(col, 0.0)})
        if not len(finite):
            notices.append(f"Sample {col} has no valid abundance; correlation unavailable.")
    correlation, counts = np.full((len(columns), len(columns)), np.nan), np.zeros((len(columns), len(columns)), int)
    for i in range(len(columns)):
        for j in range(i, len(columns)):
            valid = np.isfinite(selected[:, i]) & np.isfinite(selected[:, j])
            counts[i, j] = counts[j, i] = int(valid.sum())
            if valid.sum() >= 3:
                a, b = selected[valid, i], selected[valid, j]
                if np.std(a) > 0 and np.std(b) > 0:
                    correlation[i, j] = correlation[j, i] = float(stats.pearsonr(a, b).statistic)
    if not np.isfinite(correlation).all():
        notices.append("Some correlations are unavailable: at least 3 shared proteins and nonzero variance are required.")
    pca = _pca(selected, row_ids, config.get("pca_missing", "complete"))
    if pca["note"]:
        notices.append(pca["note"])
    if pca["policy"] == "mean":
        notices.append(f"PCA only: {pca['imputed_values']} missing values filled with protein means. Tests and correlations use unimputed data.")
    proteins = []
    kept_indices = np.flatnonzero(~excluded)
    for local_index, original in enumerate(kept_indices):
        source = frame.iloc[original]
        protein = _text(source.get(mapping.get("protein", "")))
        gene = _text(source.get(mapping.get("gene", "")))
        item = {"row_id": int(row_ids[local_index]), "protein": protein, "gene": gene}
        for group in dict.fromkeys(groups):
            values = selected[local_index, np.array(groups) == group]
            item[f"n_{group}"] = int(np.isfinite(values).sum())
            item[f"cv_{group}"] = _raw_cv(values)
        proteins.append(item)
    lookup = {item["row_id"]: item for item in proteins}
    loadings = [{"row_id": rid, "protein": lookup[rid]["protein"], "gene": lookup[rid]["gene"],
                 "PC1": loading[0], "PC2": loading[1]} for rid, loading in zip(pca["feature_row_ids"], pca["loadings"])]
    scores = [{"sample": col, "group": group, "PC1": score[0], "PC2": score[1]}
              for col, group, score in zip(columns, groups, pca["scores"])]
    summary = {"source_name": Path(path).name, "sheet": chosen, "header_row": header + 1,
               "total_rows": len(frame), "eligible_rows": total, "excluded_rows": int(excluded.sum()),
               "sample_count": len(columns), "group_counts": {g: groups.count(g) for g in dict.fromkeys(groups)},
               "abundance_scale": config.get("abundance_scale", "raw"), "normalization": config.get("normalization", "none"),
               "normalization_offsets_log2": offsets, "pca_available": pca["available"], "pca_policy": pca["policy"],
               "pca_protein_rows": pca["feature_count"], "pca_imputed_values": pca["imputed_values"],
               "explained_variance_ratio": pca["explained_variance_ratio"], "pca_note": pca["note"],
               "correlation_method": "Pearson; pairwise complete; minimum 3 proteins",
               "cv_definition": "Across-replicate sample SD / mean on linear abundance reconstructed from normalized log2 values; n >= 2."}
    summary["source_sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    summary["created_utc"] = datetime.now(timezone.utc).isoformat()
    summary["versions"] = {"numpy": np.__version__, "scipy": scipy.__version__}
    for group in dict.fromkeys(groups):
        cvs = [r[f"cv_{group}"] for r in proteins if r[f"cv_{group}"] is not None]
        summary[f"median_cv_{group}"] = float(np.median(cvs)) if cvs else None
    standard_colors = {"control": COLORS["down"], "treatment": COLORS["up"], "ungrouped": COLORS["neutral"]}
    group_colors = {group: standard_colors[group] for group in dict.fromkeys(groups) if group in standard_colors}
    remaining_colors = [color for color in CATEGORICAL if color not in group_colors.values()]
    for index, group in enumerate(group for group in dict.fromkeys(groups) if group not in group_colors):
        group_colors[group] = remaining_colors[index % len(remaining_colors)]
    plots = {
        "missingness": {"data": [{"type": "bar", "x": columns, "y": [r["missing_fraction"] for r in metrics],
                                   "customdata": [[r["missing"], r["eligible_protein_rows"]] for r in metrics],
                                   "marker": {"color": [group_colors[group] for group in groups], "line": {"color": "#111111", "width": .7}},
                                   "hovertemplate": "%{x}<br>Missing: %{customdata[0]} / %{customdata[1]}<br>%{y:.1%}<extra></extra>"}],
                        "layout": _layout("Missing values", xaxis={"title": {"text": "Sample"}}, yaxis={"title": {"text": "Missing fraction"}, "tickformat": ".0%", "range": [0, 1]}, showlegend=False)},
        "distribution": {"data": [{"type": "box", "name": col, "y": _clean(selected[:, index].tolist()),
                                     "boxpoints": False, "marker": {"color": group_colors[groups[index]]},
                                     "fillcolor": group_colors[groups[index]], "line": {"color": "#111111", "width": 1},
                                     "whiskerwidth": .5} for index, col in enumerate(columns)],
                         "layout": _layout("Abundance distribution", xaxis={"title": {"text": "Sample"}}, yaxis={"title": {"text": "log<sub>2</sub> abundance"}}, showlegend=False)},
        "correlation": {"data": [{"type": "heatmap", "x": columns, "y": columns, "z": _clean(correlation.tolist()),
                                   "customdata": counts.tolist(), "zmin": -1, "zmax": 1, "colorscale": DIVERGING, "reversescale": False,
                                   "colorbar": {"title": {"text": "Pearson r"}, "thickness": 16, "outlinecolor": "#111111", "outlinewidth": .7},
                                   "hoverongaps": False, "hovertemplate": "%{x} / %{y}<br>r = %{z:.4f}<br>Shared proteins: %{customdata}<extra></extra>"}],
                        "layout": _layout("Sample correlation", xaxis={"title": {"text": "Sample"}}, yaxis={"title": {"text": "Sample"}, "autorange": "reversed"})},
    }
    if pca["available"]:
        traces = []
        for group_index, group in enumerate(dict.fromkeys(groups)):
            group_scores = [r for r in scores if r["group"] == group]
            traces.append({"type": "scatter", "mode": "markers+text", "name": group,
                           "x": [r["PC1"] for r in group_scores], "y": [r["PC2"] for r in group_scores],
                           "text": [r["sample"] for r in group_scores], "textposition": "top center",
                           "marker": {"size": 10, "color": group_colors[group], "symbol": ["circle", "square", "diamond", "triangle-up", "triangle-down", "cross", "x"][group_index % 7], "line": {"color": "#111111", "width": .6}},
                           "hovertemplate": "%{text}<br>PC1=%{x:.3f}<br>PC2=%{y:.3f}<extra></extra>"})
        ratios = pca["explained_variance_ratio"]
        plots["pca"] = {"data": traces, "layout": _layout("PCA (centered, unscaled)",
                        xaxis={"title": {"text": f"PC1 ({ratios[0]:.1%})"}}, yaxis={"title": {"text": f"PC2 ({ratios[1]:.1%})"}})}
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    _write_table(output / "sample_metrics.csv", metrics, list(metrics[0]))
    correlation_rows = [{"sample_a": a, "sample_b": b, "pearson_r": _clean(correlation[i, j]), "n_shared": int(counts[i, j])}
                        for i, a in enumerate(columns) for j, b in enumerate(columns)]
    _write_table(output / "correlation.csv", correlation_rows, ["sample_a", "sample_b", "pearson_r", "n_shared"])
    _write_table(output / "correlation_counts.csv", correlation_rows, ["sample_a", "sample_b", "n_shared"])
    _write_table(output / "pca_loadings.csv", loadings, ["row_id", "protein", "gene", "PC1", "PC2"])
    _write_table(output / "pca_scores.csv", scores, ["sample", "group", "PC1", "PC2"])
    _write_table(output / "protein_cv.csv", proteins, ["row_id", "protein", "gene"] + [key for g in dict.fromkeys(groups) for key in (f"n_{g}", f"cv_{g}")])
    files = {"quality": "quality.html", "sample_metrics": "sample_metrics.csv", "correlation": "correlation.csv",
             "correlation_counts": "correlation_counts.csv", "pca_loadings": "pca_loadings.csv",
             "pca_scores": "pca_scores.csv", "protein_cv": "protein_cv.csv", "quality_data": "quality.json"}
    result = {"summary": summary, "warnings": notices, "files": files, "plots": plots,
              "data": {"sample_metrics": metrics, "correlation": _clean(correlation.tolist()), "correlation_counts": counts.tolist(),
                       "pca": pca, "pca_scores": scores, "protein_cv": proteins}}
    (output / "quality.json").write_text(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2), encoding="utf-8")
    _quality_html(output / "quality.html", result)
    return result


def _quality_html(path, result):
    summary = result["summary"]
    def esc(value):
        return html.escape(str(value))
    metrics = result["data"]["sample_metrics"]
    rows = "".join("<tr>" + "".join(f"<td>{esc('—' if r.get(k) is None else format(r[k], '.4g') if isinstance(r[k], float) else r[k])}</td>"
                                  for k in ("sample", "group", "observed", "missing", "missing_fraction", "median_log2", "normalization_offset_log2")) + "</tr>" for r in metrics)
    captions = {
        "missingness": "Fraction of eligible protein rows with missing abundance in each sample.",
        "distribution": "Boxes show the first quartile, median and third quartile. Whiskers extend to the most extreme values within 1.5 interquartile ranges; outlier points are not drawn.",
        "correlation": "Pairwise Pearson correlation on log2 abundance. Blank cells lack at least three shared proteins or have zero variance.",
        "pca": "Samples projected onto centered, unscaled protein abundances. Axis labels give the explained variance. Color and shape identify the sample group.",
    }
    cards = "".join(f'<section class="figure-panel"><div class="figure-actions"><strong>{chr(97 + index)}</strong>'
                    f'<button type="button" data-save-figure="{name}">Download SVG</button></div>'
                    f'<div id="{name}" class="chart"></div><p class="figure-caption">{esc(captions[name])}</p></section>'
                    for index, name in enumerate(result["plots"]))
    cv_formula = math_markup(r"\mathrm{CV}=\frac{s_A}{\overline{A}},\qquad s_A=\sqrt{\frac{\sum_{i=1}^{n}(A_i-\overline{A})^2}{n-1}},\quad n\geq 2", display=True)
    pca_formula = math_markup(r"X_c=U\Sigma V^{\mathsf{T}},\qquad \mathrm{scores}=U\Sigma", display=True)
    notices = "".join(f"<li>{esc(n)}</li>" for n in result["warnings"])
    payload = json.dumps(result["plots"], ensure_ascii=False, allow_nan=False).replace("<", "\\u003c")
    text = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Sample QC: {esc(summary['source_name'])}</title><style>
body{{margin:0;background:#fff;color:#111;font:15px/1.6 Arial,sans-serif}}body,input,button,select,textarea,pre,code,svg text{{font-family:Arial,sans-serif}}main{{max-width:1200px;margin:auto;padding:28px}}
h1{{font-size:28px;color:#111}}section,details{{background:white;border-top:1px solid #111;margin:24px 0;padding:20px 0}}summary{{cursor:pointer;color:#111}}.chart{{height:470px;min-width:0}}
table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{text-align:left;border-bottom:1px solid #aaa;padding:9px}}th{{border-bottom:2px solid #111}}.scroll{{overflow:auto}}
code{{padding:2px 5px}}a{{color:#0072B2}}.figure-actions{{display:flex;justify-content:space-between;align-items:center;gap:12px}}.figure-actions strong{{font-size:18px}}button{{font:14px Arial,sans-serif;padding:9px 13px;color:#111;border:1px solid #111;background:white;cursor:pointer}}.figure-caption{{font-size:14px;color:#111;margin:10px 0 0}}.math-display{{overflow:auto}}@media(max-width:600px){{main{{padding:18px}}.chart{{height:420px}}}}@media print{{button{{display:none}}section{{break-inside:avoid}}}}</style><main>
<h1>Sample QC</h1><p>{esc(summary['source_name'])}</p><p>{summary['sample_count']} samples, {summary['eligible_rows']} protein rows</p>
<details><summary>Methods</summary><p>Input: <code>{esc(summary['abundance_scale'])}</code>. Normalization: <code>{esc(summary['normalization'])}</code>.</p>
<p>Raw abundance ≤ 0 is missing. All values use the log2 scale. QC exclusions follow the selected flags; median offsets use included rows.</p>
<p>PCA: <code>{esc(summary['pca_policy'])}</code>, {summary['pca_protein_rows']} protein rows, centered per protein without variance scaling. Correlations use complete pairs. Differential tests do not use imputed values.</p>
{pca_formula}<p>Rows of the centered matrix correspond to samples and columns to proteins. PCA separation is not a significance test.</p>
{cv_formula}<p>CV uses reconstructed linear abundance within each group. See protein_cv.csv.</p>
<ul>{notices}</ul><p><a href="https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html">NumPy SVD</a><br><a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html">SciPy Pearson correlation</a></p></details>
{cards}<section><h2>Samples</h2><div class="scroll"><table><thead><tr><th>Sample</th><th>Group</th><th>Observed</th><th>Missing</th><th>Missing fraction</th><th>Median log2</th><th>Offset</th></tr></thead><tbody>{rows}</tbody></table></div></section></main>
<script>{get_plotlyjs()}</script><script>const figures={payload}; for(const [id,figure] of Object.entries(figures)) Plotly.newPlot(id,figure.data,figure.layout,{{responsive:true,displaylogo:false,modeBarButtonsToRemove:['sendChartToCloud','sendDataToCloud'],toImageButtonOptions:{{format:'svg',filename:id,width:1000,height:700}}}});document.querySelectorAll('[data-save-figure]').forEach(button=>button.onclick=()=>Plotly.downloadImage(document.getElementById(button.dataset.saveFigure),{{format:'svg',filename:button.dataset.saveFigure,width:1000,height:700}}));</script>{math_assets()}</html>'''
    Path(path).write_text(text, encoding="utf-8")


def analyze_enrichment(rows, gmt_text, identifier="gene", background="tested", output_dir=None, min_size=1, max_size=None):
    """One-sided hypergeometric ORA over unique, unambiguously identified genes.

    The universe is the actual valid differential-test population, not the union
    of gene sets. Duplicate rows collapse; a gene is a hit if any eligible row is
    significant. Ambiguous protein groups cannot be assigned to individual genes
    and are excluded. Identifier matching is exact and case-sensitive.
    """
    if background != "tested":
        raise ValueError("Enrichment background must be tested: valid, nonexcluded differential tests.")
    if identifier not in {"gene", "protein"}:
        raise ValueError("Use gene or protein identifiers matching the GMT namespace.")
    if isinstance(min_size, bool) or not isinstance(min_size, int) or min_size < 1:
        raise ValueError("Minimum set size must be an integer ≥ 1.")
    if max_size is not None and (isinstance(max_size, bool) or not isinstance(max_size, int) or max_size < min_size):
        raise ValueError("Maximum set size must be an integer ≥ the minimum.")
    if not isinstance(gmt_text, str) or not gmt_text.strip():
        raise ValueError("Provide a nonempty GMT gene-set file.")
    terms, warnings = {}, []
    for line_number, line in enumerate(gmt_text.splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 3 or not fields[0].strip():
            raise ValueError(f"Invalid GMT line {line_number}: use tab-separated name, description and at least one identifier.")
        name, description = fields[0].strip(), fields[1].strip()
        members = {token.strip() for token in fields[2:] if token.strip()}
        if not members:
            raise ValueError(f"GMT line {line_number} has no valid identifiers.")
        if name in terms:
            terms[name]["members"].update(members)
            warnings.append(f"Duplicate GMT set {name}: members merged and tested once.")
        else:
            terms[name] = {"description": description, "members": members}
    if not terms:
        raise ValueError("No valid gene sets found in GMT.")
    universe, hits = set(), set()
    counts = {"input_rows": len(rows), "invalid_or_excluded_rows": 0, "missing_test_rows": 0,
              "missing_identifier_rows": 0, "ambiguous_identifier_rows": 0, "eligible_rows": 0}
    for row in rows:
        if row.get("valid") is not True or row.get("excluded"):
            counts["invalid_or_excluded_rows"] += 1
            continue
        has_probability = any((value := _number(row.get(k))) is not None and 0 <= value <= 1 for k in ("p_value", "q_value", "metric_value"))
        has_score = (score := _number(row.get("significance_score"))) is not None and score >= 0
        if not has_probability and not has_score:
            counts["missing_test_rows"] += 1
            continue
        text = _text(row.get(identifier)).strip()
        # Do not expand ambiguous protein groups into multiple independent genes.
        tokens = {t.strip() for t in re.split(r"[;,|\s]+", text) if t.strip()}
        if not tokens:
            counts["missing_identifier_rows"] += 1
            continue
        if len(tokens) > 1:
            counts["ambiguous_identifier_rows"] += 1
            continue
        gene = next(iter(tokens))
        universe.add(gene)
        counts["eligible_rows"] += 1
        if row.get("significant") is True:
            hits.add(gene)
    if not universe:
        raise ValueError("No background identifiers with valid test results.")
    population, sample = len(universe), len(hits)
    table = []
    for name, term in terms.items():
        members = term["members"] & universe
        overlap = sorted(members & hits)
        size, observed = len(members), len(overlap)
        expected = sample * size / population
        eligible = size >= min_size and (max_size is None or size <= max_size)
        table.append({"term": name, "description": term["description"], "term_size_original": len(term["members"]),
                      "term_size_tested": size, "population_size": population, "hit_population_size": sample,
                      "observed_hits": observed, "expected_hits": expected,
                      "fold_enrichment": observed / expected if expected else None,
                      "p_value": float(stats.hypergeom.sf(observed - 1, population, size, sample)) if eligible else None,
                      "q_value": None, "eligible": eligible, "hit_identifiers": overlap,
                      "tested_identifiers": sorted(members)})
    qs = bh_adjust([row["p_value"] for row in table])
    for row, q in zip(table, qs):
        row["q_value"] = q
    table.sort(key=lambda row: (row["q_value"] is None, row["q_value"] or 0, row["p_value"] or 0, row["term"]))
    if counts["ambiguous_identifier_rows"]:
        warnings.append(f"Excluded {counts['ambiguous_identifier_rows']} ambiguous protein groups.")
    duplicate_count = counts["eligible_rows"] - population
    if duplicate_count:
        warnings.append(f"Collapsed {duplicate_count} duplicate identifier rows. Any eligible significant row defines a hit.")
    if not sample:
        warnings.append("No significant identifiers; all eligible gene sets have p = q = 1.")
    annotated = universe & set().union(*(term["members"] for term in terms.values()))
    summary = {**counts, "identifier": identifier, "background": background, "population_size": population,
               "hit_population_size": sample, "annotated_population_size": len(annotated),
               "unannotated_population_size": len(universe - annotated), "duplicate_rows_collapsed": duplicate_count,
               "term_count": len(table), "tested_term_count": sum(row["eligible"] for row in table),
               "min_size": min_size, "max_size": max_size,
               "method": "one-sided hypergeometric over-representation; BH across every eligible term, including zero-hit terms; size filters use tested members",
               "identifier_policy": "exact case-sensitive; ambiguous groups excluded; duplicate identifiers collapsed; any eligible significant row defines a hit"}
    summary.update(gmt_sha256=hashlib.sha256(gmt_text.encode("utf-8")).hexdigest(),
                   analysis_rows_sha256=hashlib.sha256(json.dumps(_clean(rows), ensure_ascii=False, sort_keys=True, allow_nan=False).encode("utf-8")).hexdigest(),
                   created_utc=datetime.now(timezone.utc).isoformat(), versions={"numpy": np.__version__, "scipy": scipy.__version__})
    result = {"summary": summary, "rows": table, "warnings": warnings,
              "audit": {"background_identifiers": sorted(universe), "hit_identifiers": sorted(hits),
                        "unannotated_identifiers": sorted(universe - annotated)}, "files": {}}
    if output_dir is not None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        result["files"] = {"enrichment": "enrichment.html", "enrichment_table": "enrichment.csv", "enrichment_audit": "enrichment.json"}
        _write_table(output / "enrichment.csv", table, list(table[0]))
        (output / "enrichment.json").write_text(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2), encoding="utf-8")
        _enrichment_html(output / "enrichment.html", result)
    return result


def enrich_from_rows(rows, gmt_text, output_dir, config=None):
    config = config or {}
    return analyze_enrichment(rows, gmt_text, identifier=config.get("identifier", "gene"),
                              output_dir=output_dir, min_size=config.get("min_size", 1), max_size=config.get("max_size"))


def _enrichment_html(path, result):
    summary = result["summary"]
    def esc(value):
        return html.escape(str(value))
    def number(value):
        return "—" if value is None else format(value, ".5g") if isinstance(value, float) else str(value)
    table = []
    for row in result["rows"]:
        values = [row["term"], row["description"], row["observed_hits"], row["term_size_tested"], row["expected_hits"],
                  row["fold_enrichment"], row["p_value"], row["q_value"], "Yes" if row["eligible"] else "No", "; ".join(row["hit_identifiers"])]
        table.append("<tr>" + "".join(f"<td>{esc(number(v))}</td>" for v in values) + "</tr>")
    warnings = "".join(f"<li>{esc(message)}</li>" for message in result["warnings"])
    audit = esc(json.dumps(display_provenance({"summary": summary, "audit": result["audit"]}), ensure_ascii=False, indent=2))
    enrichment_formula = math_markup(r"P(X\geq k)=\sum_{j=k}^{\min(K,n)}\frac{\binom{K}{j}\binom{N-K}{n-j}}{\binom{N}{n}}", display=True)
    bh_formula = math_markup(r"q_{(i)}=\min\left(1,\min_{j\geq i}\frac{m\,p_{(j)}}{j}\right)", display=True)
    text = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Gene-set enrichment</title><style>body{{margin:0;background:#fff;color:#111;font:15px/1.6 Arial,sans-serif}}
body,input,button,select,textarea,pre,code,svg text{{font-family:Arial,sans-serif}}h1{{color:#111;font-size:28px}}main{{max-width:1400px;margin:auto;padding:28px}}section,details{{background:white;border-top:1px solid #111;margin:24px 0;padding:20px 0}}summary{{cursor:pointer;color:#111}}
table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{text-align:left;border-bottom:1px solid #aaa;padding:9px;vertical-align:top}}
th{{background:#fff;border-bottom:2px solid #111;cursor:pointer}}.scroll{{overflow:auto}}input{{font:inherit;padding:10px;width:min(400px,100%);border:1px solid #111}}
pre{{white-space:pre-wrap;overflow-wrap:anywhere}}a{{color:#0072B2}}code{{overflow-wrap:anywhere}}@media(max-width:600px){{main{{padding:18px}}}}</style><main>
<h1>Gene-set enrichment</h1><p>{summary['population_size']} background IDs, {summary['hit_population_size']} significant IDs, {summary['tested_term_count']} tested sets</p>
<p><a href="enrichment.csv" download>All results CSV</a><br><a href="enrichment.json" download>Audit JSON</a></p>
<details><summary>Methods</summary><p>One-sided hypergeometric test. Background: valid, nonexcluded differential tests, including IDs without GMT annotations. IDs are exact, case-sensitive matches.</p>
{enrichment_formula}<p>N is the number of background IDs, K the number of tested IDs in the set, n the number of significant IDs, and k the observed overlap.</p>
{bh_formula}<p>p values are ordered from smallest to largest; m is the number of eligible gene sets.</p>
<p>Each ID counts once; any eligible significant row defines a hit. Ambiguous groups are excluded. BH includes all eligible sets, including zero-hit sets. Size is measured within the tested background: {summary['min_size']} to {esc(summary['max_size'] if summary['max_size'] is not None else 'unlimited')}.</p>
<p>Ineligible sets remain listed without p/q values. Enrichment is an association, not evidence of pathway activity or direct binding.</p><ul>{warnings}</ul>
<p><a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html">SciPy hypergeometric distribution</a></p></details>
<section><h2>Gene sets</h2><p><input id="search" type="search" placeholder="Filter set, description or hit ID"><span id="count"></span></p>
<div class="scroll"><table id="terms"><thead><tr><th>Set</th><th>Description</th><th>Hits</th><th>Background size</th><th>Expected hits</th><th>Enrichment</th><th>p</th><th>BH q</th><th>Tested</th><th>Hit IDs</th></tr></thead><tbody>{''.join(table)}</tbody></table></div></section>
<details><summary>Provenance and population audit</summary><p>GMT SHA-256: <code>{summary['gmt_sha256']}</code></p><pre>{audit}</pre></details></main>
<script>const table=document.getElementById('terms'),body=table.tBodies[0],search=document.getElementById('search');
function filter(){{let count=0;for(const row of body.rows){{row.hidden=!row.textContent.toLowerCase().includes(search.value.toLowerCase());if(!row.hidden)count++;}}document.getElementById('count').textContent=' '+count+' / '+body.rows.length+' sets';}}
search.addEventListener('input',filter);table.querySelectorAll('th').forEach((header,index)=>{{let ascending=false;header.addEventListener('click',()=>{{ascending=!ascending;const rows=Array.from(body.rows);rows.sort((a,b)=>{{const x=a.cells[index].textContent,y=b.cells[index].textContent;if(x==='—')return y==='—'?0:1;if(y==='—')return -1;const numeric=index>=2&&index<=7;return (numeric?Number(x)-Number(y):x.localeCompare(y))*(ascending?1:-1);}});rows.forEach(row=>body.appendChild(row));}});}});filter();</script>{math_assets()}</html>'''
    Path(path).write_text(text, encoding="utf-8")
