# iST Cell Type Annotation

## [`single_cell_reference/`](single_cell_reference)

Preprocessing and manual cluster annotation of the matched sc/snRNA-seq reference samples into the
final `cell_type202605` labels — see [`single_cell_reference/README.md`](single_cell_reference/README.md).

## `02_annotate_iST_cell_types.R`

Transfers the `cell_type202605` reference labels onto each Xenium/MERSCOPE sample by SingleR
(Wilcoxon DE-based label transfer against the matched reference), producing `SingleR_labels202605`.
Applied per sample and platform.

### Input

| Argument | Description |
|----------|-------------|
| `sample_id` | e.g. `MH0007` |
| `platform` | `"Xenium"` or `"MERSCOPE"` |
| `ref_type` | `"scRNA"` or `"snRNA"` — which reference folder to use (sample-specific, see table below) |
| `so.rds` | iST Seurat object from `../../preprocessing/01_iST.R` |
| `seu_<sample_id>.rds` | Matched sc/snRNA-seq reference from `single_cell_reference/03_scrna_annotate_cell_types.R`, with `cell_type202605` metadata |

### Output

`so.rds` (overwritten in place) — same object with an added `SingleR_labels202605` metadata column.

### Sample / platform / reference-type table

| Sample | Platform | ref_type |
|--------|----------|----------|
| MH0007 | Xenium | snRNA |
| MH0007 | MERSCOPE | snRNA |
| MH0026 | Xenium | snRNA |
| MH0026 | MERSCOPE | snRNA |
| ER_0114 | Xenium | scRNA |
| ER_0114 | MERSCOPE | scRNA |
| ER_0360 | Xenium | scRNA |
| ER_0360 | MERSCOPE | scRNA |
| TN_0177 | Xenium | scRNA |
| TN_0177 | MERSCOPE | scRNA |
| TN_0554 | Xenium | scRNA |

TN_0554 has no MERSCOPE row (no completed MERSCOPE preprocessing for this sample — consistent with
`../../preprocessing/README.md`).

Reference cells labelled `"mixture"` (unresolved sub-clustering in the reference annotation) are
excluded from the SingleR reference before label transfer.
