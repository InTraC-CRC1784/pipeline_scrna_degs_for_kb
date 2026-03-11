import os
import argparse
from pathlib import Path
from typing import Sequence, Dict, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from dataclasses import dataclass

@dataclass
class PreprocessConfig:
    adata_path: str
    out_dir: str
    cell_state_col: str
    condition_col: str
    sample_id_col: str
    cell_states_to_exclude: Sequence[str]
    cell_type_col: str
    cell_type_value: str
    normalize_target_sum: float = 1e4
    threshold: float = 0.0125
    gene_index: Optional[str] = None

def save_to_folder(df: pd.DataFrame, folder: str | Path, filename: str):
    """Ensures directory exists and saves DataFrame to CSV."""
    folder_path = Path(folder)
    folder_path.mkdir(parents=True, exist_ok=True)
    file_path = folder_path / filename
    df.to_csv(file_path)
    print(f"[INFO] Saved: {file_path}")

def correlation_analysis(adata: sc.AnnData, column_tosplitby: str, how: str) -> pd.DataFrame:
    """
    Groups by a metadata column and calculates sum or mean across genes.
    """
    data_source = adata.raw if adata.raw is not None else adata
    
    res = {}
    groups = adata.obs[column_tosplitby].unique()
    
    for group in groups:
        # FIX: Add .values to convert the Series to a Numpy array
        mask = (adata.obs[column_tosplitby] == group).values 
        
        sub_X = data_source.X[mask]
        
        if how == "mean":
            val = sub_X.mean(axis=0)
        else:
            val = sub_X.sum(axis=0)
            
        res[str(group)] = np.ravel(val)

    return pd.DataFrame(res, index=data_source.var_names)
    return pd.DataFrame(res, index=data_source.var_names)

def run_preprocessing(cfg: PreprocessConfig):
    print(f"[INFO] Loading AnnData from {cfg.adata_path}...")
    adata = sc.read_h5ad(cfg.adata_path)

    # --- 1. GENE INDEXING (CRITICAL: MUST BE FIRST) ---
    # We set the index before setting .raw so that .raw inherits the correct gene names.
    target_index_col = None

    if cfg.gene_index is not None:
        if cfg.gene_index in adata.var.columns:
            target_index_col = cfg.gene_index
        # Fallback: if index is numeric (0, 1, 2...), try to find a name column
        elif adata.var_names.astype(str).str.fullmatch(r"\d+").all():
            for col in ["features", "gene_ids", "symbols", "_index"]:
                if col in adata.var.columns and not adata.var[col].astype(str).str.fullmatch(r"\d+").all():
                    target_index_col = col
                    break

    if target_index_col is not None:
        print(f"[INFO] Setting gene index to column: '{target_index_col}'")
        adata.var_names = adata.var[target_index_col].astype(str)
        adata.var_names_make_unique()

    # --- 2. METADATA VALIDATION & CLEANING ---
    # Check if cell_type_col exists; if not, fill with default value
    if cfg.cell_type_col not in adata.obs.columns:
        print(f"[WARNING] '{cfg.cell_type_col}' missing. Creating with value '{cfg.cell_type_value}'.")
        adata.obs[cfg.cell_type_col] = cfg.cell_type_value

    # Standardize cell state names (remove special characters for downstream R/edgeR compatibility)
    if cfg.cell_state_col in adata.obs.columns:
        adata.obs[cfg.cell_state_col] = (
            adata.obs[cfg.cell_state_col]
            .astype(str)
            .str.replace(r"[^A-Za-z0-9_]+", "_", regex=True)
            .astype("category")
        )
    else:
        raise ValueError(f"Required column '{cfg.cell_state_col}' not found in adata.obs")

    # Filter out specific cell states
    if cfg.cell_states_to_exclude:
        initial_count = adata.n_obs
        adata = adata[~adata.obs[cfg.cell_state_col].isin(cfg.cell_states_to_exclude)].copy()
        print(f"[INFO] Excluded {initial_count - adata.n_obs} cells based on exclusion list.")

    # Fix .raw slot now that indices and cell-filtering are finalized
    adata.raw = adata

    # --- 3. PSEUDOBULK (RAW SUMS) ---
    # Creating a unique key for grouping (vectorized string concatenation is faster)
    adata.obs["pool_key"] = (
        adata.obs[cfg.condition_col].astype(str) + "__" +
        adata.obs[cfg.cell_state_col].astype(str) + "__" +
        adata.obs[cfg.sample_id_col].astype(str)
    )
    
    df_raw_sum = correlation_analysis(adata, "pool_key", "sum")
    save_to_folder(df_raw_sum, cfg.out_dir, "pseudobulk_whole.csv")

    # --- 4. NORMALIZATION & THRESHOLD FILTERING ---
    # Work on a copy to keep the main 'adata' raw for other uses
    norm_adata = adata.copy()
    sc.pp.normalize_total(norm_adata, target_sum=cfg.normalize_target_sum)
    sc.pp.log1p(norm_adata)

    # Use the same grouping logic for mean expression
    norm_adata.obs["Gene_Patient"] = adata.obs["pool_key"]
    df_norm_mean = correlation_analysis(norm_adata, "Gene_Patient", "mean")
    
    # Thresholding: set low values to 0 and drop genes that are 0 everywhere
    df_thresholded = df_norm_mean.copy()
    df_thresholded[np.abs(df_thresholded) < cfg.threshold] = 0
    df_clean = df_thresholded.loc[~(df_thresholded == 0).all(axis=1)]

    save_to_folder(df_clean, cfg.out_dir, "pseudobulk_filtering.csv")

    # --- 5. EXPORT SUMMARY TABLES ---
    # Map Cell State -> Cell Type
    translation = (
        adata.obs[[cfg.cell_state_col, cfg.cell_type_col]]
        .drop_duplicates()
        .set_index(cfg.cell_state_col)
    )
    save_to_folder(translation, cfg.out_dir, "cell_state_translation_table.csv")

    # Frequency Tables (Crosstabs)
    ct_counts = pd.crosstab(adata.obs[cfg.cell_type_col], adata.obs[cfg.sample_id_col])
    cs_counts = pd.crosstab(adata.obs[cfg.cell_state_col], adata.obs[cfg.sample_id_col])
    
    save_to_folder(ct_counts, cfg.out_dir, "cell_type_number.csv")
    save_to_folder(cs_counts, cfg.out_dir, "cell_state_number.csv")

    print(f"\n[SUCCESS] Preprocessing complete. Results are in: {cfg.out_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess Scanpy object for Pseudobulk analysis.")

    parser.add_argument("adata_path", help="Path to the input .h5ad file")
    parser.add_argument("--out-dir", default="./edgeR_results", help="Directory to save CSVs")
    parser.add_argument("--cell-state-col", required=True, help="Column name for cell states/clusters")
    parser.add_argument("--condition-col", required=True, help="Column name for experimental conditions")
    parser.add_argument("--sample-id-col", required=True, help="Column name for sample/patient IDs")
    parser.add_argument("--threshold", type=float, default=0.0125, help="Expression threshold for filtering")
    parser.add_argument("--cell-states-to-exclude", nargs="*", default=[], help="List of cell states to drop")
    parser.add_argument("--cell-type-col", required=True, help="Column name for broad cell types")
    parser.add_argument("--cell-type-val", default="Unknown", help="Default value if cell_type_col is missing")
    parser.add_argument("--gene_index", default=None, help="Optional: .var column to use as gene names")

    args = parser.parse_args()

    config = PreprocessConfig(
        adata_path=args.adata_path,
        out_dir=args.out_dir,
        cell_state_col=args.cell_state_col,
        condition_col=args.condition_col,
        sample_id_col=args.sample_id_col,
        cell_states_to_exclude=args.cell_states_to_exclude,
        cell_type_col=args.cell_type_col,
        cell_type_value=args.cell_type_val,
        threshold=args.threshold,
        gene_index=args.gene_index,
    )

    run_preprocessing(config)