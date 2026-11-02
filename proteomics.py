"""Local, auditable protein-table analysis. Excel inputs are never modified."""
from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import posixpath
import re
import warnings as pywarnings
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
import numpy as np
import pandas as pd
from scipy import stats
from lxml import etree

VERSION = "1.0.0-public-rc1"
P_FLOOR = 1e-300
DEFAULTS = {
    "sheet": None, "header": None, "mapping": {}, "mode": "existing",
    "fc_scale": "log2", "p_scale": "raw", "q_scale": "raw",
    "control": [], "treatment": [], "abundance_scale": "raw",
    "normalization": "none", "test_method": "welch", "x_base": 2, "y_base": 10,
    "fc_cutoff_mode": "log", "fc_cutoff": 1,
    "significance_metric": "q", "significance_mode": "probability",
    "significance_cutoff": 0.05, "direction": "both", "exclude_qc": False,
    "organism": 9606, "make_report": True, "make_volcano": True,
}

def _norm(value):
    return re.sub(r"[^a-z0-9\u4e00-\u9fff]+", "", str(value).lower())

def _clean(value):
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    try:
        if pd.isna(value):
            return None
    except (ValueError, TypeError):
        pass
    if isinstance(value, (str, int)):
        return value
    return str(value)

def _text(value):
    clean = _clean(value)
    return "" if clean is None else str(clean)

def _number(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, str):
        value = value.strip()
        if value.lower() in {"", "na", "nan", "n/a", "none", "null", "#n/a", "-"}:
            return None
        if re.fullmatch(r"[+-]?\d{1,3}(,\d{3})+(\.\d+)?([eE][+-]?\d+)?", value):
            value = value.replace(",", "")
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError, OverflowError):
        return None

def _detect_mapping(columns):
    cs = [str(c) for c in columns]
    ns = {c: _norm(c) for c in cs}
    def pick(exacts=(), patterns=(), reject=()):
        for token in exacts:
            for col in cs:
                if ns[col] == token and not any(re.search(x, col, re.I) for x in reject):
                    return col
        for pat in patterns:
            for col in cs:
                if re.search(pat, col, re.I) and not any(re.search(x, col, re.I) for x in reject):
                    return col
        return None
    protein = pick(("proteinids", "majorityproteinids", "proteinid", "proteinaccession", "accession", "accessions", "proteinaccessions", "uniprot", "uniprotid", "proteingroup", "proteingroups", "leadingrazorprotein", "leadingac", "proteins", "蛋白质id", "蛋白id", "蛋白编号"),
                   (r"^(?:master\s+)?protein\s+(?:accession|id)", r"uniprot", r"accession", r"protein\.ids?"), (r"count|number|description",))
    gene = pick(("genenames", "genename", "genes", "gene", "genesymbol", "基因", "基因名"), (r"gene[ ._-]?(?:name|symbol)",))
    description = pick(("proteinnames", "proteinname", "description", "proteindescription", "蛋白名称", "蛋白描述"), (r"description", r"protein[ ._-]?names?"))
    # Identification FDR is not differential-abundance significance.
    q_value = pick(("qvalue", "qvalues", "adjpval", "adjpvalue", "adjustedpvalue", "padj", "fdr", "bhqvalue"),
                   (r"(?:adjusted|adj[ ._-]?p|q[ ._-]?value|false.discovery|FDR)",),
                   (r"protein.*q.?value|global.*q.?value|pg\.q|identification|precursor",))
    p_value = pick(("pvalue", "pvalues", "pval", "p", "ttestpvalue", "studentsttestpvalue", "neglogpval", "neglogpvalue", "p值"),
                   (r"(?:^|[ _.-])p[ _.-]?val", r"(?:t[ -]?test|anova).*p[ -]?value", r"[−-]?log.*p[ -]?value"),
                   (r"adj|fdr|q[ ._-]?value|posterior|PEP",))
    fc = pick(("log2fc", "log2foldchange", "log2ratio", "logfc", "foldchange", "foldchanges", "fc", "ratio", "倍数变化"),
              (r"log2.*(?:fold|ratio|fc)", r"fold[ ._-]?change", r"t[ -]?test.*difference", r"abundance.*ratio", r"(?:^|[ ._-])logfc(?:$|[ ._(-])"),
              (r"variability|p.?value|count|cv|standard|sd",))
    if fc is None:
        ratio_candidates = [col for col in cs if re.search(r"^LFQ\s*(?:ratio\s*)?\(.*(?:control|ctrl)", col, re.I)]
        if len(ratio_candidates) == 1:
            fc = ratio_candidates[0]
    function = pick(("function", "functionalannotation", "gobiologicalprocess", "geneontology", "go", "pathway", "功能", "功能注释"),
                    (r"function|biological.?process|gene.?ontology|pathway|GO[ _.-]?(?:term|annotation)",))
    abundance = [col for col in cs if re.search(r"intensity|abundance|quantity|quantification|^lfq|^ibaq|reporter.*(?:corrected|intensity)|^TMT\d|^\d{3}[NC]?$", col, re.I)
                 and col not in {fc, p_value, q_value} and not re.search(r"ratio|variability|count|coverage|score|p.?value|q.?value|CV\b", col, re.I)]
    return {"protein": protein, "gene": gene, "description": description, "fold_change": fc,
            "p_value": p_value, "q_value": q_value, "abundance": abundance, "function": function}

def _header_score(row):
    mapping = _detect_mapping([_text(x) for x in row])
    score = sum(5 for key in ("protein", "gene", "fold_change", "p_value", "q_value") if mapping[key])
    return score + min(6, len(mapping["abundance"])) + (2 if mapping["description"] else 0)


def _xlsx_last_value_row(path, sheet_name):
    """Find data/formula extent without expanding formatting-only worksheet rows."""
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    package_ns = "http://schemas.openxmlformats.org/package/2006/relationships"
    parser = etree.XMLParser(resolve_entities=False, load_dtd=False, no_network=True)
    with zipfile.ZipFile(path) as archive:
        workbook = etree.fromstring(archive.read("xl/workbook.xml"), parser)
        sheets = workbook.findall(f"{{{main_ns}}}sheets/{{{main_ns}}}sheet")
        selected = next((item for item in sheets if item.get("name") == sheet_name), None)
        if selected is None:
            raise ValueError("Worksheet not found.")
        relation_id = selected.get(f"{{{rel_ns}}}id")
        relations = etree.fromstring(archive.read("xl/_rels/workbook.xml.rels"), parser)
        relation = next((item for item in relations.findall(f"{{{package_ns}}}Relationship")
                         if item.get("Id") == relation_id), None)
        if relation is None or relation.get("TargetMode") == "External":
            raise ValueError("Excel worksheet reference is not a local worksheet.")
        target = relation.get("Target", "")
        member = posixpath.normpath(target.lstrip("/") if target.startswith("/") else posixpath.join("xl", target))
        if not member.startswith("xl/"):
            raise ValueError("Invalid Excel worksheet reference.")
        cell_tag, row_tag = f"{{{main_ns}}}c", f"{{{main_ns}}}row"
        last_row = 0
        with archive.open(member) as stream:
            for _, element in etree.iterparse(stream, events=("end",), tag=(cell_tag, row_tag),
                                              resolve_entities=False, load_dtd=False, no_network=True):
                if element.tag == cell_tag:
                    value = element.find(f"{{{main_ns}}}v")
                    # Keep formulas even when their cached result is absent, and explicit inline strings.
                    has_value = (element.find(f"{{{main_ns}}}f") is not None
                                 or element.find(f"{{{main_ns}}}is") is not None
                                 or (value is not None and value.text is not None))
                    if has_value:
                        coordinate = re.search(r"(\d+)$", element.get("r", ""))
                        row_number = int(coordinate.group(1)) if coordinate else int(element.getparent().get("r", "0"))
                        if not 1 <= row_number <= 1048576:
                            raise ValueError("Excel cell row is outside the worksheet limits.")
                        last_row = max(last_row, row_number)
                else:
                    element.clear()
                    while element.getprevious() is not None:
                        del element.getparent()[0]
    return last_row

def _read(path, sheet=None, header=None):
    path = Path(path)
    ext = path.suffix.lower()
    if ext not in {".xlsx", ".xls", ".csv", ".tsv", ".txt"}:
        raise ValueError("Select an .xlsx, .xls, .csv or .tsv protein table.")
    if header is not None and (isinstance(header, bool) or int(header) < 0):
        raise ValueError("Header row must be at least 1.")
    if ext in {".xlsx", ".xls"}:
        try:
            with pd.ExcelFile(path) as book:
                sheets = list(book.sheet_names)
                chosen = sheet if sheet is not None else sheets[0]
                if isinstance(chosen, int):
                    chosen = sheets[chosen]
                if chosen not in sheets:
                    raise ValueError("Worksheet not found.")
                last_value_row = _xlsx_last_value_row(path, chosen) if ext == ".xlsx" else None
                if header is None:
                    peek = pd.read_excel(book, sheet_name=chosen, header=None, nrows=min(35, last_value_row) if last_value_row is not None else 35)
                    scores = [_header_score(row) for row in peek.itertuples(index=False, name=None)]
                    header = int(np.argmax(scores)) if scores else 0
                nrows = max(0, last_value_row - int(header) - 1) if last_value_row is not None else None
                frame = pd.read_excel(book, sheet_name=chosen, header=int(header), nrows=nrows)
                if nrows is not None:
                    # pandas trims trailing cached-empty formulas; their source row positions still exist.
                    frame = frame.reindex(range(nrows))
        except ImportError as exc:
            raise ValueError("Excel support requires openpyxl and xlrd. Run the dependency installer.") from exc
    else:
        sheets = [path.name]
        chosen = path.name
        encoding = None
        for candidate in ("utf-8-sig", "utf-16", "gb18030"):
            try:
                sample = path.read_text(encoding=candidate)
                encoding = candidate
                break
            except UnicodeError:
                continue
        if encoding is None:
            raise ValueError("Unsupported text encoding. Save the table as UTF-8 CSV.")
        sep = "\t" if ext in {".tsv", ".txt"} else ","
        try:
            sep = csv.Sniffer().sniff("\n".join(sample.splitlines()[:40]), delimiters=",\t;").delimiter
        except csv.Error:
            pass
        if header is None:
            peek_rows = list(csv.reader(io.StringIO(sample), delimiter=sep))[:35]
            scores = [_header_score(row) for row in peek_rows]
            header = int(np.argmax(scores)) if scores else 0
        frame = pd.read_csv(path, sep=sep, encoding=encoding, header=int(header), low_memory=False, skip_blank_lines=False)
    frame.columns = [str(col).strip() for col in frame.columns]
    seen, unique = {}, []
    for col in frame.columns:
        seen[col] = seen.get(col, 0) + 1
        unique.append(col if seen[col] == 1 else f"{col} [{seen[col]}]")
    frame.columns = unique
    return frame, sheets, str(chosen), int(header)

def _qc_columns(columns):
    return [str(col) for col in columns if _norm(col) in {"reverse", "potentialcontaminant", "contaminant", "onlyidentifiedbysite", "decoy", "isdecoy", "iscontaminant"}]

def inspect_file(path, sheet=None, header=None):
    """Read all rows; return mapping suggestions and a bounded preview."""
    frame, sheets, chosen, header = _read(path, sheet, header)
    mapping = _detect_mapping(frame.columns)
    fc_name = _norm(mapping["fold_change"] or "")
    fc_scale = "log2"
    if "log10" in fc_name:
        fc_scale = "log10"
    elif fc_name.startswith("ln"):
        fc_scale = "ln"
    elif mapping["fold_change"] and "log" not in fc_name and "difference" not in fc_name:
        fc_scale = "ratio"
    def scale(col):
        return "neglog10" if col and re.search(r"(?:[-−]\s*log|neg.*log|log10.*p|log.*p[ -]?value)", col, re.I) else "raw"
    notes = []
    if not mapping["protein"] and not mapping["gene"]:
        notes.append("Select a protein ID or gene column.")
    if mapping["fold_change"] and re.match(r"^logfc(?:$|[\s_(])", mapping["fold_change"], re.I):
        notes.append("The logFC base is unclear. Confirm the suggested log2 scale.")
    if mapping["p_value"] and _norm(mapping["p_value"]) in {"neglogpval", "neglogpvalue"}:
        notes.append("Confirm whether the transformed p-value column contains raw or adjusted p values.")
    if not mapping["fold_change"]:
        notes.append("No fold-change column detected. Select one or calculate from replicates.")
    if not mapping["p_value"] and not mapping["q_value"]:
        notes.append("No differential p/q column detected. Identification FDR is not differential significance.")
    return {"sheets": sheets, "sheet": chosen, "header": header,
            "columns": list(frame.columns), "preview": _clean(frame.head(8).to_dict(orient="records")),
            "row_count": len(frame), "mapping": mapping,
            "suggested": {"fc_scale": fc_scale, "p_scale": scale(mapping["p_value"]), "q_scale": scale(mapping["q_value"])},
            "qc_columns": _qc_columns(frame.columns), "warnings": notes}

def bh_adjust(values):
    """Benjamini–Hochberg, preserving missing indices and stable ties."""
    array = np.asarray([np.nan if value is None else value for value in values], dtype=float)
    result = np.full(len(array), np.nan)
    indexes = np.where(np.isfinite(array) & (array >= 0) & (array <= 1))[0]
    if len(indexes):
        order = indexes[np.argsort(array[indexes], kind="stable")]
        scaled = array[order] * len(order) / np.arange(1, len(order) + 1)
        result[order] = np.minimum(1, np.minimum.accumulate(scaled[::-1])[::-1])
    return [_clean(value) for value in result]

def _probability(value, scale):
    number = _number(value)
    if number is None:
        return None
    if scale == "neglog10":
        return 10.0 ** (-number) if number >= 0 else None
    return number if 0 <= number <= 1 else None

def _log_probability(value, scale):
    """Natural-log probability avoids false significance after score underflow."""
    number = _number(value)
    if number is None:
        return None
    if scale == "neglog10":
        return -number * math.log(10) if number >= 0 else None
    if not 0 <= number <= 1:
        return None
    return math.log(number) if number > 0 else -math.inf

def _bh_log(values):
    indexes = [i for i, value in enumerate(values) if value is not None and value <= 0]
    order = sorted(indexes, key=lambda i: values[i])
    adjusted = [None] * len(values)
    running = 0.0
    count = len(order)
    for position in range(count - 1, -1, -1):
        index = order[position]
        running = min(running, values[index] + math.log(count / (position + 1)))
        adjusted[index] = running
    return adjusted

def _log2_fc(value, scale):
    number = _number(value)
    if number is None:
        return None
    if scale == "ratio":
        return math.log2(number) if number > 0 else None
    if scale == "log2":
        return number
    if scale == "log10":
        return number * math.log2(10)
    return number / math.log(2)

def _near_ge(left, right):
    return left > right or math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-14)

def _near_le(left, right):
    return left < right or math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-300)

def _truthy_qc(value):
    return value is True or _text(value).strip().lower() in {"+", "true", "yes", "1", "1.0", "reverse", "contaminant", "decoy"}

def _validate_config(config, frame, detected):
    c = {**DEFAULTS, **(config or {})}
    c["mapping"] = {**detected, **(c.get("mapping") or {})}
    enums = {
        "mode": {"existing", "replicates"}, "test_method": {"welch", "student", "paired", "mannwhitney"}, "fc_scale": {"log2", "log10", "ln", "ratio"},
        "p_scale": {"raw", "neglog10"}, "q_scale": {"raw", "neglog10"},
        "abundance_scale": {"raw", "log2", "log10", "ln"}, "normalization": {"none", "median"},
        "fc_cutoff_mode": {"log", "ratio"}, "significance_metric": {"p", "q"},
        "significance_mode": {"probability", "score"}, "direction": {"both", "up", "down"},
    }
    for name, options in enums.items():
        if c[name] not in options:
            raise ValueError(f"Unsupported {name}: {c[name]}")
    for name in ("x_base", "y_base", "fc_cutoff", "significance_cutoff"):
        number = _number(c[name])
        if number is None:
            raise ValueError(f"{name} must be a finite number.")
        c[name] = number
    if c["x_base"] <= 1 or c["y_base"] <= 1:
        raise ValueError("Axis log bases must be greater than 1.")
    if c["fc_cutoff"] < (1 if c["fc_cutoff_mode"] == "ratio" else 0):
        raise ValueError("Fold-change threshold must be ≥ 1; absolute logFC threshold must be ≥ 0.")
    if c["significance_mode"] == "probability" and not 0 < c["significance_cutoff"] <= 1:
        raise ValueError("The p/q cutoff must be greater than 0 and at most 1.")
    if c["significance_mode"] == "score" and c["significance_cutoff"] < 0:
        raise ValueError("The −log significance cutoff must be ≥ 0.")
    if not c["make_report"] and not c["make_volcano"]:
        raise ValueError("Select a report, a volcano plot, or both.")
    for key, col in c["mapping"].items():
        if key == "abundance":
            if not isinstance(col, list):
                raise ValueError("Abundance columns must be a list of column names.")
            if len(col) != len(set(col)) or any(item not in frame.columns for item in col):
                raise ValueError("Abundance columns are missing or duplicated.")
        elif col and col not in frame.columns:
            raise ValueError(f"Column not found: {col}")
    if not c["mapping"].get("protein") and not c["mapping"].get("gene"):
        raise ValueError("Select a protein ID or gene column.")
    if c["mode"] == "existing":
        if not c["mapping"].get("fold_change"):
            raise ValueError("Select a fold-change column or calculate from replicates.")
        if c["significance_metric"] == "p" and not c["mapping"].get("p_value"):
            raise ValueError("Select a p-value column to use a p-value threshold.")
        if c["significance_metric"] == "q" and not (c["mapping"].get("q_value") or c["mapping"].get("p_value")):
            raise ValueError("Select a q-value column, or a p-value column for BH adjustment.")
    else:
        for name in ("control", "treatment"):
            cols = c[name]
            if not isinstance(cols, list) or len(cols) < 2 or len(cols) != len(set(cols)):
                raise ValueError("Each group needs at least 2 distinct replicate columns.")
            if any(col not in frame.columns for col in cols):
                raise ValueError("Replicate column not found.")
        if set(c["control"]) & set(c["treatment"]):
            raise ValueError("Control and treatment groups cannot share a column.")
    if c["mode"] == "replicates" and c["test_method"] == "paired" and len(c["control"]) != len(c["treatment"]):
        raise ValueError("Paired tests require equal group sizes, paired in column-list order.")
    organism = _text(c.get("organism", 9606)).strip()
    if organism and not re.fullmatch(r"\d{1,10}", organism):
        raise ValueError("Enter a numeric NCBI taxonomy ID, such as 9606 for human.")
    c["organism"] = organism
    return c

def _replicate_statistics(frame, c, qc_excluded):
    columns = c["control"] + c["treatment"]
    matrix = np.array([[np.nan if _number(value) is None else _number(value) for value in frame[col]] for col in columns], dtype=float).T
    invalid_counts = (~np.isfinite(matrix)).sum(axis=1)
    if c["abundance_scale"] == "raw":
        invalid_counts += (np.isfinite(matrix) & (matrix <= 0)).sum(axis=1)
        matrix[matrix <= 0] = np.nan
        with np.errstate(invalid="ignore", divide="ignore"):
            matrix = np.log2(matrix)
    elif c["abundance_scale"] == "log10":
        matrix *= math.log2(10)
    elif c["abundance_scale"] == "ln":
        matrix /= math.log(2)
    offsets = {}
    if c["normalization"] == "median":
        reference = matrix[~np.asarray(qc_excluded, dtype=bool)]
        with pywarnings.catch_warnings():
            pywarnings.simplefilter("ignore", RuntimeWarning)
            medians = np.nanmedian(reference, axis=0)
        if not np.isfinite(medians).all():
            raise ValueError("A sample has no valid abundance for median normalization.")
        adjustments = medians - float(np.median(medians))
        matrix -= adjustments
        offsets = {col: float(offset) for col, offset in zip(columns, adjustments)}
    n_control, results = len(c["control"]), []
    for row, invalid_count in zip(matrix, invalid_counts):
        control, treatment = row[:n_control], row[n_control:]
        from quality import run_test
        test = run_test(control, treatment, c["test_method"])
        item = {"n_control": test["n_control"], "n_treatment": test["n_treatment"],
                "missing_abundance": int(invalid_count), "mean_log2_control": test["mean_control"],
                "mean_log2_treatment": test["mean_treatment"], "log2_fc": test["mean_difference"],
                "p_value": test["p_value"], "test_statistic": test["statistic"], "statistic_note": test["note"]}
        results.append(item)
    return results, offsets

def _network_url(identifier, organism):
    return "https://string-db.org/cgi/network?" + urlencode({"identifiers": identifier, "species": organism})

def _csv_safe(value):
    if isinstance(value, (dict, list)):
        value = json.dumps(value, ensure_ascii=False)
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@", "\t", "\r")):
        return "'" + value
    return value

def _write_csv(path, rows, source_columns):
    keys = ["row_id", "protein", "gene", "description", "log2_fc", "log_fc", "p_value", "q_value", "metric_value", "significance_score", "significant", "direction", "valid", "excluded", "exclusion_reason", "qc_flags", "function", "network_url", "n_control", "n_treatment", "missing_abundance", "mean_log2_control", "mean_log2_treatment", "statistic_note"]
    keys += ["source: " + col for col in source_columns]
    with Path(path).open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            item = {key: _csv_safe(row.get(key, "")) for key in keys}
            item.update({"source: " + col: _csv_safe(row["source"].get(col)) for col in source_columns})
            writer.writerow(item)

def analyze_file(path, config, output_dir):
    """Analyze every source row and write portable, local-only artifacts."""
    frame, sheets, chosen, header = _read(path, (config or {}).get("sheet"), (config or {}).get("header"))
    c = _validate_config(config, frame, _detect_mapping(frame.columns))
    c.update(sheet=chosen, header=header)
    mapping = c["mapping"]
    qc_columns = _qc_columns(frame.columns)
    sources = _clean(frame.to_dict(orient="records"))
    qc_flags = [[col for col in qc_columns if _truthy_qc(source.get(col))] for source in sources]
    qc_excluded = [bool(c["exclude_qc"] and flags) for flags in qc_flags]
    notices, offsets = [], {}
    if qc_columns:
        notices.append("QC flags: " + ", ".join(qc_columns) + ("; flagged rows excluded." if c["exclude_qc"] else "; flagged rows retained."))
    if c["mode"] == "replicates":
        replicate_rows, offsets = _replicate_statistics(frame, c, qc_excluded)
        p_values = [row["p_value"] for row in replicate_rows]
        q_values = bh_adjust([None if excluded else p for excluded, p in zip(qc_excluded, p_values)])
        log_p_values = [_log_probability(p, "raw") for p in p_values]
        log_q_values = _bh_log([None if excluded else p for excluded, p in zip(qc_excluded, log_p_values)])
        q_source = "Benjamini–Hochberg on valid, non-QC-excluded " + c["test_method"] + " p values"
    else:
        replicate_rows = [{} for _ in sources]
        p_values = [_probability(source.get(mapping["p_value"]), c["p_scale"]) if mapping.get("p_value") else None for source in sources]
        log_p_values = [_log_probability(source.get(mapping["p_value"]), c["p_scale"]) if mapping.get("p_value") else None for source in sources]
        if mapping.get("q_value"):
            q_values = [_probability(source.get(mapping["q_value"]), c["q_scale"]) for source in sources]
            log_q_values = [_log_probability(source.get(mapping["q_value"]), c["q_scale"]) for source in sources]
            q_source = "provided column: " + mapping["q_value"]
            notices.append("Using supplied q values; the original analysis defines their adjustment method and test population.")
        else:
            q_values = bh_adjust([None if excluded else p for excluded, p in zip(qc_excluded, p_values)])
            log_q_values = _bh_log([None if excluded else p for excluded, p in zip(qc_excluded, log_p_values)])
            q_values = [math.exp(value) if value is not None else None for value in log_q_values]
            q_source = "Benjamini–Hochberg on valid, non-QC-excluded provided p values"
    if c["significance_metric"] == "p":
        notices.append("Significance uses unadjusted p values.")
    c["x_threshold"] = c["fc_cutoff"] if c["fc_cutoff_mode"] == "log" else math.log(c["fc_cutoff"], c["x_base"])
    c["y_threshold"] = c["significance_cutoff"] if c["significance_mode"] == "score" else -math.log(c["significance_cutoff"], c["y_base"])
    rows = []
    abundance_columns = list(dict.fromkeys(mapping.get("abundance", []) + (c["control"] + c["treatment"] if c["mode"] == "replicates" else [])))
    for index, (source, p, q, rep, flags, excluded) in enumerate(zip(sources, p_values, q_values, replicate_rows, qc_flags, qc_excluded)):
        protein = _text(source.get(mapping.get("protein")))
        gene = _text(source.get(mapping.get("gene")))
        log2_fc = rep.get("log2_fc") if c["mode"] == "replicates" else _log2_fc(source.get(mapping["fold_change"]), c["fc_scale"])
        x = log2_fc * math.log(2) / math.log(c["x_base"]) if log2_fc is not None else None
        if x is not None and not math.isfinite(x):
            x = None
        metric = q if c["significance_metric"] == "q" else p
        log_metric = log_q_values[index] if c["significance_metric"] == "q" else log_p_values[index]
        actual_score = -log_metric / math.log(c["y_base"]) if log_metric is not None else None
        y = min(-math.log(P_FLOOR) / math.log(c["y_base"]), actual_score) if actual_score is not None else None
        if y is not None and not math.isfinite(y):
            y = None
        reasons = []
        if not protein and not gene:
            reasons.append("Missing protein ID and gene")
        if x is None:
            reasons.append("Missing or invalid fold change")
        if metric is None:
            reasons.append("Missing or invalid " + c["significance_metric"] + " value")
        if rep.get("statistic_note"):
            reasons.append(rep["statistic_note"])
        if excluded:
            reasons.append("Excluded by QC flags: " + ", ".join(flags))
        valid = x is not None and y is not None and bool(protein or gene) and not excluded
        passes_p = False
        if metric is not None:
            if c["significance_mode"] == "probability":
                passes_p = log_metric is not None and _near_le(log_metric, math.log(c["significance_cutoff"]))
            else:
                passes_p = actual_score is not None and _near_ge(actual_score, c["significance_cutoff"])
        passes_direction = x is not None and (c["direction"] == "both" or (c["direction"] == "up" and x > 0) or (c["direction"] == "down" and x < 0))
        significant = bool(valid and _near_ge(abs(x), c["x_threshold"]) and passes_p and passes_direction)
        direction = ("up" if x > 0 else "down" if x < 0 else "unchanged") if x is not None else "unknown"
        identifier = (gene or protein).split(";")[0].strip()
        rows.append(_clean({"row_id": index + header + 2, "protein": protein, "gene": gene,
            "description": _text(source.get(mapping.get("description"))),
            "log2_fc": log2_fc, "log_fc": x, "x": x, "y": y, "p_value": p, "q_value": q,
            "metric_value": metric, "significance_score": actual_score, "significant": significant,
            "direction": direction, "valid": valid, "excluded": excluded, "exclusion_reason": "；".join(reasons),
            "qc_flags": flags, "abundance": {col: _number(source.get(col)) for col in abundance_columns},
            "function": _text(source.get(mapping.get("function"))),
            "network_url": _network_url(identifier, c["organism"]), "source": source, **rep}))
    valid_rows = [row for row in rows if row["valid"]]
    hits = [row for row in rows if row["significant"]]
    zero_count = sum(row["metric_value"] == 0 for row in valid_rows)
    tiny_count = sum(row["metric_value"] is not None and 0 < row["metric_value"] < P_FLOOR for row in valid_rows)
    if zero_count or tiny_count:
        notices.append(f"{zero_count} zero and {tiny_count} sub-{P_FLOOR:g} significance values. Plot floor: {P_FLOOR:g}; decisions use original values.")
    invalid_count = sum(not row["valid"] and not row["excluded"] for row in rows)
    if invalid_count:
        notices.append(f"{invalid_count} rows lack a valid identifier, fold change or significance value. Retained in CSV, omitted from the plot.")
    ids = [row["protein"] or row["gene"] for row in rows if row["protein"] or row["gene"]]
    duplicate_count = len(ids) - len(set(ids))
    if duplicate_count:
        notices.append(f"{duplicate_count} duplicate identifier rows retained without merging.")
    if c["mode"] == "replicates":
        notices.append(c["test_method"] + " test on log2 abundance; no imputation. FC is mean treatment minus mean control on the log2 scale. Paired tests and FC use complete pairs only.")
    summary = {"total_rows": len(rows), "valid_rows": len(valid_rows), "invalid_rows": invalid_count,
        "excluded_rows": sum(row["excluded"] for row in rows), "significant_rows": len(hits),
        "up": sum(row["direction"] == "up" for row in hits), "down": sum(row["direction"] == "down" for row in hits),
        "unchanged_significant": sum(row["direction"] == "unchanged" for row in hits),
        "duplicate_identifier_rows": duplicate_count, "zero_metric_rows": zero_count,
        "bh_test_count": sum(p is not None and not exc for p, exc in zip(p_values, qc_excluded)),
        "source_name": Path(path).name, "sheet": chosen, "header_row": header + 1,
        "q_source": q_source, "normalization_offsets_log2": offsets}
    manifest = {"software": "Devil's Proteomics", "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(), "source_name": Path(path).name,
        "source_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "configuration": c, "summary": summary, "warnings": notices,
        "method": {"contrast": "treatment / control", "p_test": "two-sided " + c["test_method"] + " test on log2 abundances" if c["mode"] == "replicates" else "provided results",
            "q_value": q_source, "plot_probability_floor": P_FLOOR,
            "threshold_comparison": "inclusive, floating tolerance 1e-12 relative",
            "network_links": "No online network query is performed", "local_processing": True},
        "versions": {"numpy": np.__version__, "pandas": pd.__version__}}
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {"all_rows": "all_proteins.csv", "top_hits": "top_hits.csv", "manifest": "manifest.json"}
    _write_csv(output_dir / files["all_rows"], rows, list(frame.columns))
    _write_csv(output_dir / files["top_hits"], hits, list(frame.columns))
    (output_dir / files["manifest"]).write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    result = {"summary": summary, "rows": rows, "files": files, "config": c, "warnings": notices}
    from matlab_export import write_volcano_bundle
    files["matlab_volcano"] = "volcano_matlab.zip"
    write_volcano_bundle(output_dir / files["matlab_volcano"], {**result, "manifest": manifest})
    from proteomics_reports import write_report, write_volcano
    if c["make_volcano"]:
        files["volcano"] = "volcano.html"
        write_volcano(output_dir / files["volcano"], result)
    if c["make_report"]:
        files["report"] = "report.html"
        write_report(output_dir / files["report"], result, manifest)
    return _clean(result)
