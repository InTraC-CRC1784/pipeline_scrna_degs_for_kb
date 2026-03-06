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
conda create -n pipeline_scrna_degs_for_kb -c conda-forge scanpy python-igraph leidenalg bioconda::bioconductor-edger r-dplyr
conda env export --name pipeline_scrna_degs_for_kb --from-history > environment.yaml
```
 
## Execution
```
python pipeline_main.py \
    dataset.h5ad \
    --cell-state-col cell_type_combined \
    --condition-col Conditions \
    --sample-id-col Sample_ID \
    --comparison-normal-value SR \
    --out-dir pipeline_AF_LAA_SR_cAF \
    --region LAA \
    --threshold 0.0125 \
    --cell-type-col celltype \
    --cellstates-excluded "doublets,Unknown" \
    --cell-type-val heart \
    --cell-level celltypelevel \
    --species human \
    --year 2026 \
    --paper AF_LAA_SR_cAF \
    --output-file-name edgeR_results.csv \
    --min-cells 3 \
    --min-cells-per-state 3 \
    --fdr-threshold 0.05
```

**Parameters**
- adata_filepath – Path to the input .h5ad file (required)
- out_dir – Output directory for results (required)
- cell_state_col – Column defining cell states  (required)
- condition_col – Column defining experimental conditions  (required)
- sample_id_col – Column containing sample IDs  (required)
- region – Tissue/region identifier (default Unknown)
- threshold – threshold for highly expressed genes (default 0.0125)
- cell_type_col – Column defining cell types (if not present in the data, the column will be created)
- cellstates_excluded – Cell states to exclude (Default: empty list)
- cell_type_val – (default Unknown)
- cell_level – (default Unknown)
- comparison_normal_value – Reference condition (e.g., control)  (required)
- species – Species name (default Unknown)
- year – Study year (default Unknown)
- paper – Study identifier (default Unknown)
- output_file_name – Name of result file (default results.csv)
- min_cells – Minimum cells per sample (for Wilcoxon, default 3)
- min_cells_per_state – Minimum cells per cell state (for edgeR, default 3)
- fdr-treshold - cutoff for filtering

## Output
The output is a csv file in canonical form for ingestion with the InTraC knowledge graph adapter.

## Example

Using the sepsis data as an example.
```
python pipeline_main.py ../data/kaiser_sepsis/sepsis_final_with_celltypes_raw.h5ad --cell-state-col cell_type --condition-col condition --sample-id-col hash.ID --comparison-normal-value control --out-dir ../data/kaiser_sepsis/deg_pipeline/ --region blood --cell-type-col cell_type --species Human --year 2024 --paper Kaiser_SciAdv_2024 --output-file-name edger_results.csv --fdr-threshold 0.05
```
