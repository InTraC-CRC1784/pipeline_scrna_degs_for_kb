import os
import subprocess
import scanpy as sc

def _add_arg(cmd, key, value):
    """
    Add '--key=value' only if value exists
    """
    if value is None:
        return
    if isinstance(value, str) and value.strip() == "":
        return
    cmd.append(f"--{key}={value}")

def run_preprocessing(
        adata_filepath,
        cell_state_col,
        condition_col,
        sample_id_col,
        out_dir="pipeline_out",
        threshold=0.0125,
        cellstates_excluded="",
        cell_type_col="celltype",
        cell_type_val="Unknown",
        gene_index=None
):
    preprocess_cmd = [
        "python",
        "01_count_data_preprocessing_edgeR.py",
        adata_filepath,
    ]
    _add_arg(preprocess_cmd, "out-dir", out_dir)
    _add_arg(preprocess_cmd, "cell-state-col", cell_state_col)
    _add_arg(preprocess_cmd, "condition-col", condition_col)
    _add_arg(preprocess_cmd, "sample-id-col", sample_id_col)
    _add_arg(preprocess_cmd, "threshold", threshold)
    _add_arg(preprocess_cmd, "cell-states-to-exclude", cellstates_excluded)
    _add_arg(preprocess_cmd, "cell-type-col", cell_type_col)
    _add_arg(preprocess_cmd, "cell-type-val", cell_type_val)
    _add_arg(preprocess_cmd, "gene_index", gene_index)

    subprocess.run(preprocess_cmd, check=True)

def run_edgeR(
    cell_state_col,
    condition_col,
    sample_information_col,
    comparison_normal_value,
    pb_whole_path,
    pb_filtering_path,
    cell_state_translation_path,
    cell_state_number_path,
    cell_type_col="celltype",
    cell_level="Unknown",
    region="Unknown",
    species="Unknown",
    year="Unknown",
    study_name="Unknown",    # Changed from paper
    study_id="Unknown",      
    disease="Unknown",        
    disease_id="Unknown",    
    tissue="Unknown",         
    tissue_id="Unknown",      
    ontology_path=None,       
    gene_mapping_path=None,   
    results_dir="pipeline_out",
    threshold=0.0125,
    min_cells_per_state=3,
    output_file_name="edgeR.csv",
    fdr_threshold=1
):
    r_cmd = [
        "Rscript",
        "02_manual_test_edgeR.r",
    ]

    _add_arg(r_cmd, "Cell_Type_Col", cell_type_col)
    _add_arg(r_cmd, "Cell_State_Col", cell_state_col)
    _add_arg(r_cmd, "CELL_LEVEL", cell_level)
    _add_arg(r_cmd, "REGION", region)
    _add_arg(r_cmd, "Condition_Col", condition_col)
    _add_arg(r_cmd, "Sample_Information_Col", sample_information_col)
    _add_arg(r_cmd, "comparison_normal_value", comparison_normal_value)
    
    _add_arg(r_cmd, "species", species)
    _add_arg(r_cmd, "study_year", str(year))
    _add_arg(r_cmd, "study_name", study_name) # Passing study_name to R script
    _add_arg(r_cmd, "study_id", study_id)
    _add_arg(r_cmd, "disease", disease)
    _add_arg(r_cmd, "disease_id", disease_id)
    _add_arg(r_cmd, "tissue", tissue)
    _add_arg(r_cmd, "tissue_id", tissue_id)
    
    _add_arg(r_cmd, "ontology_path", ontology_path)
    _add_arg(r_cmd, "gene_mapping_path", gene_mapping_path)

    _add_arg(r_cmd, "RESULTS_DIR", results_dir)
    _add_arg(r_cmd, "pb_whole_path", pb_whole_path)
    _add_arg(r_cmd, "pb_filtering_path", pb_filtering_path)
    _add_arg(r_cmd, "cell_state_translation_path", cell_state_translation_path)
    _add_arg(r_cmd, "cell_state_number_path", cell_state_number_path)
    _add_arg(r_cmd, "min_cells_per_state", min_cells_per_state)
    _add_arg(r_cmd, "threshold", threshold)
    _add_arg(r_cmd, "output_file_name", output_file_name)
    _add_arg(r_cmd, "fdr_threshold", fdr_threshold)

    subprocess.run(r_cmd, check=True)

def run_full_pipeline(
    adata_filepath,
    cell_state_col,
    condition_col,
    sample_id_col,
    comparison_normal_value,
    out_dir="pipeline_out",
    region="Unknown",
    threshold=0.0125,
    cell_type_col="celltype",
    cellstates_excluded="",
    cell_type_val="Unknown",
    cell_level="Unknown",
    species="Unknown",
    year="Unknown",
    study_name="Unknown",    # Updated parameter name
    study_id="Unknown",     
    disease="Unknown",      
    disease_id="Unknown",   
    tissue="Unknown",        
    tissue_id="Unknown",    
    ontology_path=None,     
    gene_mapping_path=None, 
    output_file_name="results.csv",
    min_cells=3,
    min_cells_per_state=3,
    fdr_threshold=1,
    gene_index=None
):
    print("INFO: Deciding test method")
    adata = sc.read_h5ad(adata_filepath)
    sample_counts = adata.obs.groupby(condition_col)[sample_id_col].nunique()
    min_reps = sample_counts.min()

    if min_reps < 3:
        print("INFO: Min Samples < 3, start Wilcoxon")
        # Add your Wilcoxon logic here if needed
        pass 
    else:
        print("INFO: Min Samples > 2. Starting preprocessing.")
        run_preprocessing(
            adata_filepath=adata_filepath,
            out_dir=out_dir,
            cell_state_col=cell_state_col,
            condition_col=condition_col,
            sample_id_col=sample_id_col,
            threshold=threshold,
            cellstates_excluded=cellstates_excluded,
            cell_type_col=cell_type_col,
            cell_type_val=cell_type_val,
            gene_index=gene_index
        )

        pb_whole_path = os.path.join(out_dir, "pseudobulk_whole.csv")
        pb_filtering_path = os.path.join(out_dir, "pseudobulk_filtering.csv")
        cell_state_translation_path = os.path.join(out_dir, "cell_state_translation_table.csv")
        cell_state_number_path = os.path.join(out_dir, "cell_state_number.csv")

        print("INFO: Starting edgeR")
        run_edgeR(
            cell_type_col=cell_type_col,
            cell_state_col=cell_state_col,
            cell_level=cell_level,
            region=region,
            condition_col=condition_col,
            sample_information_col=sample_id_col,
            comparison_normal_value=comparison_normal_value,
            species=species,
            year=year,
            study_name=study_name, # Passed correctly now
            study_id=study_id,
            disease=disease,
            disease_id=disease_id,
            tissue=tissue,
            tissue_id=tissue_id,
            ontology_path=ontology_path,
            gene_mapping_path=gene_mapping_path,
            results_dir=out_dir,
            pb_whole_path=pb_whole_path,
            pb_filtering_path=pb_filtering_path,
            cell_state_translation_path=cell_state_translation_path,
            cell_state_number_path=cell_state_number_path,
            threshold=threshold,
            output_file_name=output_file_name,
            min_cells_per_state=min_cells_per_state,
            fdr_threshold=fdr_threshold
        )
        print("INFO: edgeR finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run full scRNAseq pipeline")
    parser.add_argument("adata_filepath")
    parser.add_argument("--cell-state-col", required=True)
    parser.add_argument("--condition-col", required=True)
    parser.add_argument("--sample-id-col", required=True)
    parser.add_argument("--comparison-normal-value", required=True)
    parser.add_argument("--out-dir", default="pipeline_out")
    parser.add_argument("--region", default="Unknown")
    parser.add_argument("--threshold", type=float, default=0.0125)
    parser.add_argument("--cell-type-col", default="celltype")
    parser.add_argument("--cellstates-excluded", default="")
    parser.add_argument("--cell-type-val", default="Unknown")
    parser.add_argument("--cell-level", default="Unknown")
    parser.add_argument("--species", default="Human")
    parser.add_argument("--year", default="Unknown")
    parser.add_argument("--study-name", default="Unknown") # Updated CLI flag
    
    parser.add_argument("--study-id", default="Unknown")
    parser.add_argument("--disease", default="Unknown")
    parser.add_argument("--disease-id", default="Unknown")
    parser.add_argument("--tissue", default="Unknown")
    parser.add_argument("--tissue-id", default="Unknown")
    parser.add_argument("--ontology-path", default=None)
    parser.add_argument("--gene-mapping-path", default=None)

    parser.add_argument("--output-file-name", default="results.csv")
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--min-cells-per-state", type=int, default=3)
    parser.add_argument("--fdr-threshold", type=float, default=1)
    parser.add_argument("--gene-index", default=None)

    args = parser.parse_args()

    run_full_pipeline(
        adata_filepath=args.adata_filepath,
        cell_state_col=args.cell_state_col,
        condition_col=args.condition_col,
        sample_id_col=args.sample_id_col,
        comparison_normal_value=args.comparison_normal_value,
        out_dir=args.out_dir,
        region=args.region,
        threshold=args.threshold,
        cell_type_col=args.cell_type_col,
        cellstates_excluded=args.cellstates_excluded,
        cell_type_val=args.cell_type_val,
        cell_level=args.cell_level,
        species=args.species,
        year=args.year,
        study_name=args.study_name, # Updated argument name
        study_id=args.study_id,
        disease=args.disease,
        disease_id=args.disease_id,
        tissue=args.tissue,
        tissue_id=args.tissue_id,
        ontology_path=args.ontology_path,
        gene_mapping_path=args.gene_mapping_path,
        output_file_name=args.output_file_name,
        min_cells=args.min_cells,
        min_cells_per_state=args.min_cells_per_state,
        fdr_threshold=args.fdr_threshold,
        gene_index=args.gene_index
    )