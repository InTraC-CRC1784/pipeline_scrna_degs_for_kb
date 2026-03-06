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

def run_preprocessing(adata_filepath, 
                      cell_state_col, 
                      condition_col, 
                      sample_id_col, 
                      out_dir = "pipeline_out",    #set default 
                      threshold = 0.0125, 
                      cellstates_excluded = "", 
                      cell_type_col = "celltype", 
                      cell_type_val = "Unknown" ):
    preprocess_cmd = [
        "python",
        "01_count_data_preprocessing_edgeR.py",
        adata_filepath,
    ]
    
    #add arguments to arg parser
    _add_arg(preprocess_cmd, "out-dir", out_dir)
    _add_arg(preprocess_cmd, "cell-state-col", cell_state_col)
    _add_arg(preprocess_cmd, "condition-col", condition_col)
    _add_arg(preprocess_cmd, "sample-id-col", sample_id_col)
    _add_arg(preprocess_cmd, "threshold", threshold)
    _add_arg(preprocess_cmd, "cell-states-to-exclude", cellstates_excluded)
    _add_arg(preprocess_cmd, "cell-type-col", cell_type_col)
    _add_arg(preprocess_cmd, "cell-type-val", cell_type_val)

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
    cell_type_col = "celltype",
    cell_level = "Unknown",
    region = "Unknown",
    species = "Unknown",
    year = "Unknown",
    paper = "Unknown",
    results_dir = "pipeline_out",    
    threshold = 0.0125,
    min_cells_per_state = 3,
    output_file_name  = "edgeR.csv",
    fdr_threshold = 1
):
    
    rscript_exe="Rscript"
    r_cmd = [
        rscript_exe,
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
    _add_arg(r_cmd, "year", str(year))
    _add_arg(r_cmd, "paper", paper)
    _add_arg(r_cmd, "RESULTS_DIR", results_dir)

    _add_arg(r_cmd, "pb_whole_path", pb_whole_path)
    _add_arg(r_cmd, "pb_filtering_path", pb_filtering_path)
    _add_arg(r_cmd, "cell_state_translation_path", cell_state_translation_path)
    _add_arg(r_cmd, "cell_state_number_path", cell_state_number_path)
    _add_arg(r_cmd, "min_cells_per_state", min_cells_per_state)
    _add_arg(r_cmd, "threshold", str(threshold))
    _add_arg(r_cmd, "output_file_name", output_file_name)
    _add_arg(r_cmd, "fdr_threshold", fdr_threshold)

    subprocess.run(r_cmd, check=True)
    
def run_wilcoxon(adata_path,                  
                 cell_state_col,
                 condition_col,
                 sample_id_col,
                 reference_value,
                 out_dir = "pipeline_out",
                 cell_type_col = "celltype",
                 cell_type_val = "Unknown",
                 min_cells = 3,
                 region = "Unknown",
                 annotation_level = "Unknown",
                 species = "Unknown",
                 paper = "Unknown",
                 year = "Unknown",
                 fdr_threshold = 1
                ):
    wilcoxon_cmd = [
        "python",
        "wilcoxon.py",
        adata_path,
    ]

    _add_arg(wilcoxon_cmd, "out-dir", out_dir)
    _add_arg(wilcoxon_cmd, "cell-state-col", cell_state_col)
    _add_arg(wilcoxon_cmd, "condition-col", condition_col)
    _add_arg(wilcoxon_cmd, "sample-id-col", sample_id_col)
    _add_arg(wilcoxon_cmd, "reference-value", reference_value)
    _add_arg(wilcoxon_cmd, "cell-type-col", cell_type_col)
    _add_arg(wilcoxon_cmd, "cell-type-val", cell_type_val)
    _add_arg(wilcoxon_cmd, "min-cells", min_cells)
    _add_arg(wilcoxon_cmd, "region", region)
    _add_arg(wilcoxon_cmd, "annotation-level", annotation_level)
    _add_arg(wilcoxon_cmd, "species", species)
    _add_arg(wilcoxon_cmd, "paper", paper)
    _add_arg(wilcoxon_cmd, "year", str(year))
    _add_arg(wilcoxon_cmd, "fdr_threshold", fdr_threshold)

    subprocess.run(wilcoxon_cmd, check = True)




def run_full_pipeline(
    adata_filepath,    
    cell_state_col,
    condition_col,
    sample_id_col,
    comparison_normal_value,
    out_dir = "pipeline_out",
    region = "Unknown",
    threshold = 0.0125,
    cell_type_col = "celltype",
    cellstates_excluded = "",
    cell_type_val = "Unknown",
    cell_level = "Unknown",
    species = "Unknown",
    year = "Unknown",
    paper = "Unknown",
    output_file_name = "results.csv",
    min_cells = 3,
    min_cells_per_state = 3,
    fdr_threshold = 1
):
    """
    Decision between Wilcoxon (<3 samples per condition) 
    or preprocessing & edgeR

    """
    print("INFO: Deciding test method")

    adata = sc.read_h5ad(adata_filepath)

    sample_counts = (
        adata.obs
        .groupby(condition_col)[sample_id_col]
        .nunique()
    )


    min_reps = sample_counts.min()

    if min_reps < 3:
        print("INFO: Min Samples < 3, start Wilcoxon")

        run_wilcoxon(
            adata_path=adata_filepath,
            out_dir=out_dir,
            cell_state_col=cell_state_col,
            condition_col=condition_col,
            sample_id_col=sample_id_col,
            reference_value=comparison_normal_value,
            cell_type_col=cell_type_col,
            cell_type_val=cell_type_val,
            region=region,
            annotation_level=cell_level,
            species=species,
            paper=paper,
            year=year,
            min_cells = min_cells,
            fdr_threshold=fdr_threshold
        )

        print("INFO: Wilcoxon finished")
        return

    else:
        print("INFO: Min Samples > 2")

        print("INFO: start preprocessing")
        run_preprocessing(
            adata_filepath=adata_filepath,
            out_dir=out_dir,
            cell_state_col=cell_state_col,
            condition_col=condition_col,
            sample_id_col=sample_id_col,
            threshold=str(threshold),
            cellstates_excluded=cellstates_excluded,
            cell_type_col=cell_type_col,
            cell_type_val=cell_type_val
        )

        # get pathes for edgeR
        pb_whole_path = os.path.join(out_dir, "pseudobulk_whole.csv")
        pb_filtering_path = os.path.join(out_dir, "pseudobulk_filtering.csv")
        cell_state_translation_path = os.path.join(out_dir, "cell_state_translation_table.csv")
        cell_state_number_path = os.path.join(out_dir, "cell_state_number.csv")
        results_dir = out_dir

        print("INFO: start edgeR")
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
            paper=paper,
            results_dir=results_dir,
            pb_whole_path=pb_whole_path,
            pb_filtering_path=pb_filtering_path,
            cell_state_translation_path=cell_state_translation_path,
            cell_state_number_path=cell_state_number_path,
            threshold=threshold,
            output_file_name=output_file_name,
            min_cells_per_state = min_cells_per_state,
            fdr_threshold=fdr_threshold
        )

        print("INFO: edgeR finished")
        return


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run full scRNAseq pipeline")

    # required arguments
    parser.add_argument("adata_filepath", type=str)
    parser.add_argument("--cell-state-col", required=True)
    parser.add_argument("--condition-col", required=True)
    parser.add_argument("--sample-id-col", required=True)
    parser.add_argument("--comparison-normal-value", required=True)

    # optional
    parser.add_argument("--out-dir", default="pipeline_out")
    parser.add_argument("--region", default="Unknown")
    parser.add_argument("--threshold", type=float, default=0.0125)
    parser.add_argument("--cell-type-col", default="celltype")
    parser.add_argument("--cellstates-excluded", default="")
    parser.add_argument("--cell-type-val", default="Unknown")
    parser.add_argument("--cell-level", default="Unknown")
    parser.add_argument("--species", default="Unknown")
    parser.add_argument("--year", default="Unknown")
    parser.add_argument("--paper", default="Unknown")
    parser.add_argument("--output-file-name", default="results.csv")
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--min-cells-per-state", type=int, default=3)
    parser.add_argument("--fdr-threshold", type=float, default=1)

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
        paper=args.paper,
        output_file_name=args.output_file_name,
        min_cells=args.min_cells,
        min_cells_per_state = args.min_cells_per_state,
        fdr_threshold = args.fdr_threshold
    )

