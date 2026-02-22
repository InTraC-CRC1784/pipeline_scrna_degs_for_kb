import os
from pathlib import Path
from typing import Optional, Sequence, Dict

import numpy as np
import pandas as pd
import scanpy as sc

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

@dataclass()
class PreprocessConfig:
    adata_path: str
    out_dir: str
    cell_state_col: str
    condition_col: str
    sample_id_col: str
    cell_states_to_exclude: Sequence[str]
    cell_type_col: str
    cell_type_value: str
    normalize_target_sum:  1e4 
    threshold: float 


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)

def paste0(i, j, k):
    return(str(i)+"__"+str(j)+"__"+str(k))

def save_to_folder(df, folder, filename):
    # Check if the folder exists, and create it if it doesn't
    if not os.path.exists(folder):
        os.makedirs(folder)
    file_path = os.path.join(folder, filename)
    df.to_csv(file_path)


def correlation_analysis(scanpy_object, column_tosplitby, sum_or_mean):
    d = {}
    for cluster_number in np.unique(scanpy_object.obs[column_tosplitby].values):
        scanpy_object_subset = scanpy_object[scanpy_object.obs[column_tosplitby].isin([cluster_number])]

        if sum_or_mean == "mean":
            d[cluster_number] = np.squeeze(np.asarray(scanpy_object_subset.raw.X.mean(axis=0)))  # np.log(+1)
        elif sum_or_mean == "sum":
            d[cluster_number] = np.squeeze(np.asarray(scanpy_object_subset.raw.X.sum(axis=0)))  # np.log(+1)
        del scanpy_object_subset

    return d

def _pseudobulk_by_key(adata: sc.AnnData, key_col: str, how: str) -> pd.DataFrame:

    if adata.raw is None:
        raise ValueError("adata.raw is None. Set adata.raw = adata before calling.")

    groups = pd.unique(adata.obs[key_col].values)
    d: Dict[str, np.ndarray] = {}
    for g in groups:
        sub = adata[adata.obs[key_col].isin([g])]
        X = sub.raw.X
        if how == "sum":
            vec = np.squeeze(np.asarray(X.sum(axis=0)))
        elif how == "mean":
            vec = np.squeeze(np.asarray(X.mean(axis=0)))
        else:
            raise ValueError("how must be 'sum' or 'mean'")
        d[str(g)] = vec

    df = pd.DataFrame(d, index=adata.raw.var_names)
    return df

def run_preprocessing(cfg: PreprocessConfig) -> dict[str, pd.DataFrame]:

    adata_raw = sc.read_h5ad(str(cfg.adata_path)) 

    if adata_raw.var_names.astype(str).str.fullmatch(r"\d+").all():
        for col in ["features", "_index"]:
            if col in adata_raw.var.columns:
                # Only replace if the column actually contains non-numeric gene names
                if not adata_raw.var[col].astype(str).str.fullmatch(r"\d+").all():
                    adata_raw.var_names = adata_raw.var[col].astype(str)
                    adata_raw.var_names_make_unique()
                    print(f"[INFO] Fixed numeric var_names using column '{col}'")
                    break

    # Standardize cell state names
    if cfg.cell_state_col in adata_raw.obs.columns:         
        adata_raw.obs[cfg.cell_state_col] = (
            adata_raw.obs[cfg.cell_state_col]
            .astype(str)
            .str.replace(r"[^A-Za-z0-9_]+", "_", regex=True)  
            .astype("category")
        )

    # Optional cell state exclusions
    if cfg.cell_states_to_exclude and cfg.cell_state_col in adata_raw.obs.columns:
        mask = ~adata_raw.obs[cfg.cell_state_col].isin(list(cfg.cell_states_to_exclude))
        adata_raw = adata_raw[mask].copy()


    if cfg.cell_type_col in adata_raw.obs.columns:
        pass   
    else: adata_raw.obs[cfg.cell_type_col] = cfg.cell_type_value  


    global_all = adata_raw
    global_all.raw = global_all


    sample_col_pool = cfg.sample_id_col

    sample_col_gp = cfg.sample_id_col
    if sample_col_gp not in global_all.obs.columns and "Sample_ID" in global_all.obs.columns:
        sample_col_gp = "Sample_ID"


    global_all.obs['pool_key'] = [paste0(i, j, k) for i, j, k in zip(global_all.obs[cfg.condition_col],
                                                                     global_all.obs[cfg.cell_state_col],
                                                                     global_all.obs[cfg.sample_id_col]
                                                                     )]

    global_subset = global_all.X[0:100, 0:100]

    global_all.raw = global_all
    x = correlation_analysis(global_all, 'pool_key', 'sum')
    x = pd.DataFrame(x)
    x.index = global_all.var.index

    save_to_folder(x, cfg.out_dir, 'pseudobulk_whole.csv')
    global_all.obs['Gene_Patient']=[str(i) + "__" + str(j) + "__" + str(k) for i, j, k in zip(
                                                            global_all.obs[cfg.condition_col],  # 42
                                                               global_all.obs[cfg.cell_state_col],
                                                               global_all.obs[sample_col_gp])]

    subset_ = sc.AnnData(
        global_all[:, global_all.var.index].raw.X,  #
        obs=global_all.obs,
        var=global_all[:, global_all.var.index].raw.var
    )

    sc.pp.normalize_total(subset_, target_sum=cfg.normalize_target_sum)
    sc.pp.log1p(subset_)
    
    subset_.raw = subset_
    x = correlation_analysis(subset_, 'Gene_Patient', 'mean') 
    x = pd.DataFrame(x)
    x.index = global_all.var.index

    save_to_folder(x, cfg.out_dir, 'pseudobulk_filtering.csv')

    df_thresholded = x.copy()

    df_thresholded[np.abs(df_thresholded) < cfg.threshold] = 0

    df_clean = df_thresholded.loc[~(df_thresholded == 0).all(axis=1)]

    pseudobulk_filtering = _pseudobulk_by_key(subset_, "Gene_Patient", how="mean")
    save_to_folder(pseudobulk_filtering, Path(cfg.out_dir), "pseudobulk_filtering.csv")

    x = pd.Series(global_all.obs[cfg.cell_type_col].values,
                  index=global_all.obs[cfg.cell_state_col]).to_dict()
    x = pd.DataFrame.from_dict(x, orient='index')

    save_to_folder(x, cfg.out_dir, 'cell_state_translation_table.csv')

    pd.crosstab(global_all.obs[cfg.cell_type_col],
                global_all.obs[cfg.cell_state_col]
                )

    x = pd.crosstab(global_all.obs[cfg.cell_type_col],
                     global_all.obs[cfg.sample_id_col])

    y = pd.crosstab(global_all.obs[cfg.cell_state_col],
                    global_all.obs[cfg.sample_id_col])

    save_to_folder(x, cfg.out_dir, 'cell_type_number.csv')
    save_to_folder(y, cfg.out_dir, 'cell_state_number.csv')



if __name__ == "__main__":

    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("adata_path", help="Path to .h5ad")
    p.add_argument("--out-dir", default="./edgeR_results")
    p.add_argument("--cell-state-col")
    p.add_argument("--condition-col")
    p.add_argument("--sample-id-col")
    p.add_argument("--threshold", type=float, default=0.0125)
    p.add_argument("--cell-states-to-exclude", nargs="*",default=[])
    p.add_argument("--cell-type-col", type = str)
    p.add_argument("--cell-type-val", type = str)
    args = p.parse_args()

    cfg = PreprocessConfig(
        adata_path=args.adata_path,
        out_dir=args.out_dir,

        cell_state_col=args.cell_state_col,
        condition_col=args.condition_col,
        sample_id_col=args.sample_id_col,

        cell_states_to_exclude=args.cell_states_to_exclude,

        cell_type_col=args.cell_type_col,
        cell_type_value=args.cell_type_val,

        normalize_target_sum=1e4,
        threshold=args.threshold,
    )

    run_preprocessing(cfg)