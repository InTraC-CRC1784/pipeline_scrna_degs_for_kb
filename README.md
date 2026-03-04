This code creates a canonical DE table which can be used as input for a knowledgegraph. 
It can be called like this:

pipeline_main.run_full_pipeline(
    adata_filepath='data.h5ad',
    out_dir='./pipeline_out',
    cell_state_col='celltype',
    condition_col='Conditions',
    sample_id_col='Sample_ID',
    region='LAA',
    threshold=0.0125,
    cell_type_col='celltype',
    cellstates_excluded='',
    cell_type_val='Unknown',
    cell_level='Unknown',
    comparison_normal_value='control',
    species='human',
    year='2026',
    paper='A_FIB',
    output_file_name='results.csv',
    min_cells=3,
    min_cells_per_state=3
)

Parameters

adata_filepath – Path to the input .h5ad file (required)

out_dir – Output directory for results (required)

cell_state_col – Column defining cell states  (required)

condition_col – Column defining experimental conditions  (required)

sample_id_col – Column containing sample IDs  (required)

region – Tissue/region identifier (default Unknown)

threshold – threshold for highly expressed genes (default 0.0125)

cell_type_col – Column defining cell types (required)

cellstates_excluded – Cell states to exclude (Default: empty list)

cell_type_val – (default Unknown)

cell_level – (default Unknown)

comparison_normal_value – Reference condition (e.g., control)  (required)

species – Species name (default Unknown)

year – Study year (default Unknown)

paper – Study identifier (default Unknown)

output_file_name – Name of result file (default results.csv)

min_cells – Minimum cells per sample (for Wilcoxon, default 3)

min_cells_per_state – Minimum cells per cell state (for edgeR, default 3)
