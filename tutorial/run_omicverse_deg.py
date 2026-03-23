#!/usr/bin/env python3
"""
Run OmicVerse differential expression analysis and export DENetwork-ready files.

Example
-------
python tutorial/run_omicverse_deg.py \
    --config tutorial/configs/influenza_article_standard.json

or

python tutorial/run_omicverse_deg.py \
    --input DENetwork/example/GSE192528_RawCountsAnnotated.xlsx \
    --name influenza \
    --output-dir tutorial/output_influenza \
    --gene-column Gene \
    --treatment-pattern influenza \
    --control-pattern uninfected \
    --treatment-label influenza \
    --control-label uninfected
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument(
        "--config",
        default=None,
        help="Path to a JSON config file. CLI arguments override config values.",
    )
    config_args, _ = config_parser.parse_known_args()

    config = {}
    config_path = None
    if config_args.config:
        config_path = Path(config_args.config).expanduser().resolve()
        with config_path.open("r", encoding="utf-8") as handle:
            config = json.load(handle)
        if not isinstance(config, dict):
            raise SystemExit("Config file must contain a single JSON object.")

        allowed_keys = {
            "input",
            "name",
            "output_dir",
            "gene_column",
            "sheet_name",
            "treatment_pattern",
            "control_pattern",
            "treatment_label",
            "control_label",
            "n_cpus",
            "padj_threshold",
            "fc_threshold",
            "fold_threshold",
            "logp_max",
            "min_base_mean",
            "min_log2_base_mean",
        }
        unknown_keys = sorted(set(config) - allowed_keys)
        if unknown_keys:
            raise SystemExit(
                "Unknown config keys: " + ", ".join(unknown_keys)
            )

    parser = argparse.ArgumentParser(
        description="Run OmicVerse bulk differential expression analysis.",
    )
    parser.add_argument(
        "--config",
        default=config_args.config,
        help="Path to a JSON config file. CLI arguments override config values.",
    )
    parser.add_argument(
        "--input",
        required="input" not in config,
        help="Input count matrix file (.xlsx, .tsv, or .csv).",
    )
    parser.add_argument(
        "--name",
        required="name" not in config,
        help="Dataset name used in output filenames, for example 'influenza'.",
    )
    parser.add_argument(
        "--output-dir",
        required="output_dir" not in config,
        help="Directory for output files.",
    )
    parser.add_argument(
        "--gene-column",
        default="Gene",
        help="Column containing gene symbols or Ensembl gene IDs. Default: Gene",
    )
    parser.add_argument(
        "--sheet-name",
        default=0,
        help="Excel sheet name or index when --input is .xlsx. Default: 0",
    )
    parser.add_argument(
        "--treatment-pattern",
        required="treatment_pattern" not in config,
        help="Substring used to identify treatment sample columns.",
    )
    parser.add_argument(
        "--control-pattern",
        required="control_pattern" not in config,
        help="Substring used to identify control sample columns.",
    )
    parser.add_argument(
        "--treatment-label",
        default=None,
        help="Condition label written to metadata. Default: --treatment-pattern",
    )
    parser.add_argument(
        "--control-label",
        default=None,
        help="Condition label written to metadata. Default: --control-pattern",
    )
    parser.add_argument(
        "--n-cpus",
        type=int,
        default=8,
        help="Number of CPUs passed to PyDESeq2. Default: 8",
    )
    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for DE genes. Default: 0.05",
    )
    parser.add_argument(
        "--fc-threshold",
        type=float,
        default=-1.0,
        help="Absolute log2FC threshold. Use -1 for OmicVerse auto-threshold. Default: -1",
    )
    parser.add_argument(
        "--fold-threshold",
        type=int,
        default=0,
        help="Histogram bin index used by OmicVerse auto-thresholding. Default: 0",
    )
    parser.add_argument(
        "--logp-max",
        type=int,
        default=6,
        help="Maximum plotted -log10(qvalue) used by OmicVerse. Default: 6",
    )
    parser.add_argument(
        "--min-base-mean",
        type=float,
        default=-1.0,
        help="Filter result rows with baseMean <= this value before fold-change filtering. Set to a negative value to disable. Default: -1",
    )
    parser.add_argument(
        "--min-log2-base-mean",
        type=float,
        default=1.0,
        help="Filter result rows with log2(BaseMean) <= this value before fold-change filtering. Set to a negative value to disable. Default: 1.0",
    )
    parser.set_defaults(**config)
    args = parser.parse_args()
    args._config_path = config_path
    return args


def load_omicverse():
    try:
        import omicverse as ov
        _ = ov.bulk.pyDEG
    except Exception as exc:  # pragma: no cover - environment-dependent
        raise SystemExit(
            "Failed to import omicverse bulk DE tools. "
            "Make sure the omicverse environment also has pydeseq2, cytoolz, and frozendict installed.\n"
            f"Original error: {exc}"
        ) from exc

    try:
        import pydeseq2  # noqa: F401
    except Exception as exc:  # pragma: no cover - environment-dependent
        raise SystemExit(
            "PyDESeq2 is required by ov.bulk.pyDEG(method='DEseq2'). "
            "Install it first, for example with `conda install -n omicverse pydeseq2 -c conda-forge -c bioconda`.\n"
            f"Original error: {exc}"
        ) from exc

    return ov


def read_input_table(path: Path, sheet_name: str | int, gene_column: str) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        table = pd.read_excel(path, sheet_name=sheet_name)
    elif suffix in {".tsv", ".txt"}:
        table = pd.read_csv(path, sep="\t")
    elif suffix == ".csv":
        table = pd.read_csv(path)
    else:
        raise SystemExit(f"Unsupported input format: {path.suffix}")

    if gene_column not in table.columns:
        raise SystemExit(f"Gene column '{gene_column}' was not found in {path}")

    return table


def find_sample_columns(
    table: pd.DataFrame,
    gene_column: str,
    treatment_pattern: str,
    control_pattern: str,
) -> tuple[list[str], list[str]]:
    treatment_columns = [
        col
        for col in table.columns
        if col != gene_column and treatment_pattern in str(col)
    ]
    control_columns = [
        col
        for col in table.columns
        if col != gene_column and control_pattern in str(col)
    ]

    overlap = sorted(set(treatment_columns) & set(control_columns))
    if overlap:
        raise SystemExit(
            "Some sample columns matched both treatment and control patterns: "
            + ", ".join(overlap)
        )
    if not treatment_columns:
        raise SystemExit(
            f"No sample columns matched treatment pattern '{treatment_pattern}'."
        )
    if not control_columns:
        raise SystemExit(
            f"No sample columns matched control pattern '{control_pattern}'."
        )

    return treatment_columns, control_columns


def build_count_matrix(
    table: pd.DataFrame,
    gene_column: str,
    sample_columns: list[str],
) -> pd.DataFrame:
    counts = table[[gene_column] + sample_columns].copy()
    counts = counts.dropna(subset=[gene_column])
    counts[gene_column] = counts[gene_column].astype(str).str.strip()
    counts = counts[counts[gene_column] != ""]

    for column in sample_columns:
        counts[column] = pd.to_numeric(counts[column], errors="coerce")

    if counts[sample_columns].isna().any().any():
        bad_columns = [
            col for col in sample_columns if counts[col].isna().any()
        ]
        raise SystemExit(
            "Found non-numeric values in count columns: " + ", ".join(bad_columns)
        )

    rounded = counts[sample_columns].round()
    if not (rounded == counts[sample_columns]).all().all():
        raise SystemExit(
            "Count matrix contains non-integer values. "
            "OmicVerse DESeq2 expects raw integer counts."
        )

    counts[sample_columns] = rounded.astype(int)
    return counts


ENSEMBL_GENE_PATTERNS = {
    "human": re.compile(r"^ENSG\d+(?:\.\d+)?$"),
    "mouse": re.compile(r"^ENSMUSG\d+(?:\.\d+)?$"),
}

GENE_MAPPING_FILES = {
    "human": "pair_GRCh38.tsv",
    "mouse": "pair_GRCm39.tsv",
}


def detect_ensembl_species(gene_values: pd.Series) -> str | None:
    genes = gene_values.astype(str).str.strip()
    genes = genes[genes != ""]
    if genes.empty:
        return None

    sample = genes.head(min(len(genes), 500))
    for species, pattern in ENSEMBL_GENE_PATTERNS.items():
        if sample.str.match(pattern).mean() >= 0.8:
            return species
    return None


@contextlib.contextmanager
def pushd(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def ensure_gene_mapping_file(ov, output_dir: Path, species: str) -> Path:
    mapping_filename = GENE_MAPPING_FILES[species]
    mapping_path = output_dir / "genesets" / mapping_filename
    if mapping_path.exists():
        return mapping_path

    with pushd(output_dir):
        ov.utils.download_geneid_annotation_pair()

    if not mapping_path.exists():
        raise SystemExit(
            f"Expected OmicVerse mapping file was not created: {mapping_path}"
        )
    return mapping_path


def maybe_map_gene_ids(
    ov,
    count_matrix: pd.DataFrame,
    output_dir: Path,
) -> tuple[pd.DataFrame, str | None, Path | None]:
    detected_species = detect_ensembl_species(pd.Series(count_matrix.index))
    if detected_species is None:
        return count_matrix, None, None

    mapping_path = ensure_gene_mapping_file(ov, output_dir, detected_species)

    mapped = count_matrix.copy()
    mapped.index = (
        pd.Index(mapped.index)
        .astype(str)
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
    )
    mapped = ov.bulk.Matrix_ID_mapping(
        mapped,
        str(mapping_path),
        keep_unmapped=False,
    )
    return mapped, detected_species, mapping_path


def write_core_result_table(result: pd.DataFrame, output_path: Path) -> None:
    core = result[
        ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    ].copy()
    core.index.name = "Gene"
    core.to_csv(output_path, sep="\t")


def write_de_gene_table(result: pd.DataFrame, output_path: Path, sig_label: str | None) -> None:
    if sig_label is None:
        subset = result.loc[result["sig"] != "normal"].copy()
    else:
        subset = result.loc[result["sig"] == sig_label].copy()

    subset = subset[
        ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    ].copy()
    subset.index.name = "Gene"
    subset.to_csv(output_path, sep="\t")


def main() -> int:
    args = parse_args()
    ov = load_omicverse()

    config_dir = args._config_path.parent if args._config_path else None

    def resolve_path(value: str) -> Path:
        path = Path(value).expanduser()
        if not path.is_absolute() and config_dir is not None:
            path = (config_dir / path).resolve()
        else:
            path = path.resolve()
        return path

    input_path = resolve_path(args.input)
    output_dir = resolve_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    treatment_label = args.treatment_label or args.treatment_pattern
    control_label = args.control_label or args.control_pattern

    table = read_input_table(input_path, args.sheet_name, args.gene_column)
    treatment_columns, control_columns = find_sample_columns(
        table=table,
        gene_column=args.gene_column,
        treatment_pattern=args.treatment_pattern,
        control_pattern=args.control_pattern,
    )
    sample_columns = control_columns + treatment_columns
    count_table = build_count_matrix(
        table=table,
        gene_column=args.gene_column,
        sample_columns=sample_columns,
    )
    raw_count_rows = len(count_table)

    metadata = pd.DataFrame(
        {
            "sample": sample_columns,
            "condition": [control_label] * len(control_columns)
            + [treatment_label] * len(treatment_columns),
        }
    )
    metadata.to_csv(output_dir / f"{args.name}_meta_data.tsv", sep="\t", index=False)

    count_matrix = count_table.set_index(args.gene_column)
    count_matrix, detected_ensembl_species, mapping_path = maybe_map_gene_ids(
        ov=ov,
        count_matrix=count_matrix,
        output_dir=output_dir,
    )

    count_output = output_dir / f"{args.name}_count_data.tsv"
    count_matrix.reset_index().rename(
        columns={"index": args.gene_column}
    ).to_csv(count_output, sep="\t", index=False)

    dds = ov.bulk.pyDEG(count_matrix)
    deduplicated_counts = dds.drop_duplicates_index()

    normalized = ov.bulk.deseq2_normalize(deduplicated_counts)
    normalized.index.name = "Gene"
    normalized.to_csv(
        output_dir / f"{args.name}_normalized_matrix.tsv",
        sep="\t",
    )

    result = dds.deg_analysis(
        treatment_columns,
        control_columns,
        method="DEseq2",
        n_cpus=args.n_cpus,
    ).copy()

    if args.min_base_mean >= 0:
        result = result.loc[result["baseMean"] > args.min_base_mean].copy()

    if args.min_log2_base_mean >= 0:
        result = result.loc[result["log2(BaseMean)"] > args.min_log2_base_mean].copy()

    dds.result = result

    dds.foldchange_set(
        fc_threshold=args.fc_threshold,
        pval_threshold=args.padj_threshold,
        logp_max=args.logp_max,
        fold_threshold=args.fold_threshold,
    )
    result = dds.result.copy()

    write_core_result_table(
        result,
        output_dir / f"res_matrix_{args.name}.tsv",
    )

    gene_log2fc = result["abs(log2FC)"].copy()
    gene_log2fc.to_csv(
        output_dir / f"gene_log2fc_{args.name}.tsv",
        sep="\t",
        header=False,
    )

    write_de_gene_table(
        result,
        output_dir / f"DE_all_{args.name}.tsv",
        sig_label=None,
    )
    write_de_gene_table(
        result,
        output_dir / f"DE_pos_{args.name}.tsv",
        sig_label="up",
    )
    write_de_gene_table(
        result,
        output_dir / f"DE_neg_{args.name}.tsv",
        sig_label="down",
    )

    print("Input:", input_path)
    print("Output directory:", output_dir)
    print("Treatment samples:", treatment_columns)
    print("Control samples:", control_columns)
    print("Raw count rows:", raw_count_rows)
    if detected_ensembl_species is None:
        print("Gene ID mapping:", "not required")
    else:
        print(
            "Gene ID mapping:",
            f"detected {detected_ensembl_species} Ensembl IDs and mapped them with {mapping_path}",
        )
        print("Rows after gene ID mapping:", len(count_matrix))
    print("Deduplicated rows used for DE:", len(deduplicated_counts))
    print("Rows in res_matrix:", len(result))
    print("Minimum baseMean filter:", args.min_base_mean)
    print("Minimum log2(BaseMean) filter:", args.min_log2_base_mean)
    print("Auto/selected log2FC threshold:", getattr(dds, "fc_max", None))
    print("DE_all:", int((result["sig"] != "normal").sum()))
    print("DE_pos:", int((result["sig"] == "up").sum()))
    print("DE_neg:", int((result["sig"] == "down").sum()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
