# Unified pipeline to ingest cell type specific differential expression into the InTraC knowledge base

The pipeline starts with annotated single cell data in anndata objects.
Depending on the number of samples present in the file the pipeline performs either pseudo bulk differential gene expression (DEG) analysis using edgeR (for more than 3 samples per condition) or using the Wilcoxon rank rum test (if less than 3 samples per condition are present).

## Installation
The pipeline depends on the following packages:
- python
- scanpy
- R
- edgeR

You can install the dependencies using conda
```
conda create -n pipeline_scrna_degs_for_kb environment.yaml
conda activate pipeline_scrna_degs_for_kb
```

To create the config file use these commands
```
conda create -n pipeline_scrna_degs_for_kb -c conda-forge scanpy python-igraph leidenalg bioconda::bioconductor-edger r-dplyr r-tidyr
conda env export --name pipeline_scrna_degs_for_kb --from-history | grep -v "prefix:" > environment.yaml
```
 
## Execution
```
python your_script_name.py \
    --adata-filepath "dataset.h5ad" \
    --cell-state-col "cell_type" \
    --condition-col "condition" \
    --sample-id-col "sample_ID" \
    --comparison-normal-value "control" \
    --out-dir "results_edgeR" \
    --region "blood" \
    --threshold 0.0125 \
    --cell-type-col "cell_origin" \
    --cell-type-val "blood" \
    --species "Human" \
    --year "2024" \
    --study-name "NAME" \
    --study-id "author_disease" \
    --disease "disease" \
    --disease-id "MONDO:" \
    --tissue "tissue" \
    --tissue-id "UBERON:" \
    --ontology-path "ontology_mapping.csv" \
    --gene-mapping-path "gene_symbol_to_id.csv" \
    --output-file-name "edgeR_results_final.csv"
```

**Parameters**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `adata_filepath` | Path to the input `.h5ad` file | — |
| `out_dir` | Output directory where results will be saved | — |
| `cell_state_col` | Column name defining cell states | — |
| `condition_col` | Column name defining experimental conditions | — |
| `sample_id_col` | Column name containing sample identifiers | — |
| `comparison_normal_value` | Reference condition used for comparisons (e.g., `control`) | — |
| `region` | Tissue or anatomical region identifier | `Unknown` |
| `threshold` | Threshold used to define highly expressed genes | `0.0125` |
| `cell_type_col` | Column defining cell types. If not present in the dataset, the column will be created | — |
| `cell_type_val` | Default value assigned to the cell type column if created | `Unknown` |
| `cellstates_excluded` | List of cell states to exclude from the analysis | `[]` |
| `cell_level` | Cell-level annotation or hierarchy level | `Unknown` |
| `species` | Species name | `Unknown` |
| `year` | Year of the study | `Unknown` |
| `study_name` | Study identifier or name | `Unknown` |
| `output_file_name` | Name of the output results file | `results.csv` |
| `min_cells` | Minimum number of cells per sample required for the Wilcoxon test | `3` |
| `min_cells_per_state` | Minimum number of cells per cell state required for edgeR analysis | `3` |
| `fdr_threshold` | False discovery rate cutoff used for filtering results | — |
| `disease` | Name of the disease being studied | — |
| `disease_id` | MONDO identifier for the disease | — |
| `tissue` | Tissue from which samples were collected | — |
| `tissue_id` | UBERON identifier for the tissue | — |
| `ontology_path` | Path to a CSV file containing cell ontology mappings | — |
| `gene_mapping_path` | Path to a file mapping gene names to Ensembl IDs | — |

## Output
The output is a csv file in canonical form for ingestion with the InTraC knowledge graph adapter.

## Example

Using the sepsis data as an example.

Download the data from NEFELI using this link:

Execute using the yaml metadata file from NEFELI (replace the file path for the pipeline accordingly):

```
python ~/work/intrac_agent/pipeline_scrna_degs_for_kb/pipeline_main.py --metadata-file kaiser_sepsis.yaml
```

Or execute using command line parameters:

```
DATASET_PATH=../data/kaiser_sepsis/
python pipeline_main.py --adata-filepath $DATASET_PATH/sepsis_final_with_celltypes_raw.h5ad \ 
  --cell-state-col cell_type --condition-col condition --sample-id-col hash.ID \ 
  --comparison-normal-value control --out-dir $DATASET_PATH/deg_pipeline/ \ 
  --region blood --cell-type-col cell_type --species Human --year 2024 \ 
  --paper Kaiser_SciAdv_2024 --output-file-name edger_results.csv --fdr-threshold 0.05
```
