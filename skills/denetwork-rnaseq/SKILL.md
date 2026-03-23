---
name: denetwork-rnaseq
description: Prepare RNA-seq count tables for OmicVerse, choose upregulated or downregulated differential-expression targets, determine human or mouse priors, curate receptor lists when missing, and run OmicVerse plus DENetwork end to end. Use when a user wants to go from a raw count matrix or transcriptomics table to DENetwork results with minimal manual setup.
---

# DENetwork RNA-seq Workflow

Run this workflow when the user wants a full RNA-seq to DENetwork pipeline, not just a one-off command.

## Working assumptions

- Prefer the existing DENetwork repository in the workspace. Locate it by finding `DENetwork/run_denetwork.py`.
- Prefer the existing `omicverse` conda environment. If it is missing, state that clearly and create or repair it before running anything else.
- Treat the user's input matrix as raw counts unless the file clearly contains normalized or log-transformed values.
- Ask concise direct questions only when the answer cannot be inferred safely.

## Required user inputs

Do not start OmicVerse until these items are known or safely inferred:

1. Input count matrix.
2. Control samples and treatment samples.
3. Differential-expression target set: `up` or `down`.
4. Species: `human` or `mouse`.
5. Receptor strategy: provided receptor file or automatic curation.

Treat sample grouping as the only truly mandatory biological input. Species and receptor strategy can often be inferred or auto-generated, but control versus treatment must be correct.

## Question order

When information is missing, ask in this order:

1. `Which samples are controls and which are treatment or experiment?`
2. `Use upregulated targets or downregulated targets?`
3. `Is this dataset human or mouse?`
4. `Do you already have a receptor list, or should I generate one from CellTalkDB?`

## Infer the count-table format

- Open the file and inspect the first few rows and column names.
- Accept `.tsv`, `.txt`, `.csv`, `.xlsx`, or `.xls`.
- Require one gene-symbol column plus one column per sample.
- Require integer counts for the sample columns.
- If the table is already in the format expected by `tutorial/run_omicverse_deg.py`, use it directly.
- If the table is close but not quite correct, reshape it into a wide count matrix with:
  - one gene column
  - one sample column per replicate
  - integer count values

If the matrix is not raw counts, stop and tell the user that OmicVerse DESeq2 expects raw integer counts.

## Infer the sample groups

- Inspect sample names for common group tokens such as `WT`, `KO`, `Ctrl`, `Control`, `Case`, `Treat`, `Treated`, `Disease`, `Normal`, `Tumor`, or explicit condition names.
- If exactly two groups are obvious, use them.
- If the direction of comparison matters, treat the perturbation or disease group as treatment and the baseline group as control.
- If the grouping is ambiguous, stop and ask one concise question naming the candidate columns or patterns you found.
- Never guess sample groups when the column names do not support a clear mapping.

## Choose the DE target direction

If the user did not specify whether they care about upregulated or downregulated DE genes, ask:

`Use upregulated targets or downregulated targets?`

Map the answer to DENetwork inputs like this:

- `up` -> use `DE_pos_<name>.tsv`
- `down` -> use `DE_neg_<name>.tsv`

Remember that current DENetwork node weights use `|log2FC|`, so direction is controlled by which DE file you pass as targets, not by signed node weights.

## Choose the species

If the user did not specify species, infer it from gene symbols and context:

- Prefer `mouse` when many genes look like mouse symbols such as `Alox15`, `Pcna`, `Fth1`, `0610007P14Rik`, or other mixed-case MGI-style names.
- Prefer `human` when most symbols are HGNC-style uppercase.
- Use filenames, dataset descriptions, or nearby notes if they make the species obvious.

If inference is still uncertain, ask:

`Is this dataset human or mouse?`

Use the species to select the DENetwork prior:

- `human` -> default prior `DENetwork/data/ppi_ptm_pd_hgnc.txt`
- `mouse` -> run `DENetwork/run_denetwork.py --species mouse`, which uses `DENetwork/data/ppi_biogrid_mouse.txt`

## Choose the receptor strategy

Use this priority order:

1. A receptor file explicitly provided by the user.
2. An obvious receptor file already present in the repo for this dataset.
3. Automatic CellTalkDB receptor curation.

If the user already provided a receptor file, do not replace it with CellTalkDB.

If no receptor list is available, auto-curate it:

- Download the species-specific CellTalkDB ligand-receptor table.
- Extract the `receptor_gene_symbol` column.
- Uppercase the symbols to match DENetwork's internal normalization.
- Intersect them with the genes present in the expression matrix.
- Write the result to a temporary receptor file and use that file for DENetwork.

Use these CellTalkDB downloads:

- `mouse_lr_pair.txt`
- `human_lr_pair.txt`

The download endpoint is the CellTalkDB `handler/download.php` form target, not a stable static file URL. Recreate the POST request if needed.

## Run OmicVerse differential expression

Prefer the repository script:

- `tutorial/run_omicverse_deg.py`

Use article-style defaults unless the user explicitly requests different thresholds:

- `baseMean > 50`
- `padj < 0.05`
- `|log2FC| >= 0.6`

Build the command from the inferred file format and sample patterns. Typical outputs go into:

- `.tmp/<dataset>_omicverse/`

Confirm that the script produces at least:

- `DE_pos_<name>.tsv`
- `DE_neg_<name>.tsv`
- `DE_all_<name>.tsv`
- `gene_log2fc_<name>.tsv`

## Run DENetwork

Prefer the repository script:

- `DENetwork/run_denetwork.py`

Use:

- the selected DE target file based on `up` or `down`
- the generated `gene_log2fc_<name>.tsv`
- the receptor file
- `--species mouse` for mouse data
- a sensible worker count such as `-P 8` or `-P 16` depending on the machine

Use a model name that captures dataset, direction, and species, for example:

- `alox15_mouse_down`
- `influenza_human_up`

## Validate the run

After OmicVerse:

- report the inferred treatment and control groups
- report the number of DE targets in `DE_pos`, `DE_neg`, and `DE_all`

After DENetwork:

- report the size of the initialized graph
- report the number of receptors and targets that actually entered the graph
- report the final result paths:
  - `DENetwork/results/<model_name>/optimal_network.sif`
  - `DENetwork/results/<model_name>/nodes_ranking.tsv`
  - `DENetwork/results/<model_name>/genes_go.txt`

If the run fails, inspect the log first and fix the concrete blocker instead of stopping at a generic error message.

## Default execution pattern

When the user asks to run a new dataset, follow this order:

1. Inspect the input table.
2. Confirm or infer control and treatment samples.
3. Confirm or ask about DE direction.
4. Confirm or infer species.
5. Use a provided receptor list or auto-curate one from CellTalkDB.
6. Run OmicVerse.
7. Run DENetwork.
8. Summarize the outputs and key counts.
