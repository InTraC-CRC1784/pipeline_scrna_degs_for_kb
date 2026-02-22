from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, Dict
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import harmonypy as hm

@dataclass()
class WilcoxonConfig():
    data_path: str
    out_dir: str

    # columns in adata.obs
    cell_state_col: str
    condition_col: str
    sample_id_col: str

    comparison_normal_value : str
    cell_states_to_exclude: Sequence[str]

    cell_type_col: str
    cell_type_value: str  
    threshold: float  
    min_cells: int
    
    region: str
    annotation_level: str
    species: str
    paper: str
    year: str

def _run_one_wilcoxon(
    adata_sub,
    condition_col: str,
    case_value: str,
    reference_value: str) -> pd.DataFrame:

    print(f"[INFO]: running wilcoxon for case: '{case_value}', reference: '{reference_value}'")

    if adata_sub.raw is None:
        adata_sub.raw = adata_sub

    sc.tl.rank_genes_groups(
        adata_sub,
        groupby=condition_col,
        groups=[case_value],
        reference=reference_value,
        method="wilcoxon",
        use_raw=True,
    )

    df = sc.get.rank_genes_groups_df(adata_sub, group=case_value)
    return df

def run_Wilcoxon(cfg: WilcoxonConfig) -> Dict[str, pd.DataFrame]:
    out_base = Path(cfg.out_dir)
    out_base.mkdir(parents=True, exist_ok=True)

    all_results = []
    adata = sc.read_h5ad(cfg.data_path)

    if cfg.cell_type_col and cfg.cell_type_col not in adata.obs.columns:
        adata.obs[cfg.cell_type_col] = cfg.cell_type_value

    exclude = set(cfg.cell_states_to_exclude or [])
    cell_states = [
        cs for cs in adata.obs[cfg.cell_state_col].astype(str).unique().tolist()
        if cs not in exclude
    ]

    out_base = Path(cfg.out_dir)
    results: Dict[str, pd.DataFrame] = {}

    for cs in cell_states:
        ad_cs = adata[adata.obs[cfg.cell_state_col].astype(str) == str(cs)].copy()
        conds = ad_cs.obs[cfg.condition_col].astype(str).unique().tolist()

        if cfg.comparison_normal_value not in conds:
            continue

        other_times = [c for c in conds if c != cfg.comparison_normal_value]
        dfs = []

        for case in other_times:
            counts = ad_cs.obs[cfg.condition_col].astype(str).value_counts()
            n_case = int(counts.get(str(case), 0))
            n_ref = int(counts.get(str(cfg.comparison_normal_value), 0))

            if n_case < cfg.min_cells or n_ref < cfg.min_cells:
                print(
                    f"[SKIP] cell_state={cs} case={case} ref={cfg.comparison_normal_value} "
                    f"(n_case={n_case}, n_ref={n_ref}, min_cells={cfg.min_cells})"
                )
                continue
                
            df = _run_one_wilcoxon(
                adata_sub=ad_cs,
                condition_col=cfg.condition_col,
                case_value=case,
                reference_value=cfg.comparison_normal_value,
            )

            
            df["cell_type"] = str(cs)  
            df["case"] = str(case)
            df["reference"] = str(cfg.comparison_normal_value)
            df["cell_type_value"] = cfg.cell_type_value
            df["comparison"] = df["reference"] + "_" + df["case"]    
            df["test"] = "wilcoxon"       
            df["Region"] = cfg.region
            df["annotation_level"] = cfg.annotation_level
            df["species"] = cfg.species
            df["paper"] = cfg.paper
            df["year"] = cfg.year
     
            df = df.rename(columns={
                "names": "Gene",
                "logfoldchanges": "logFC",
                "pvals": "PValue",
                "pvals_adj": "FDR_all"})
            
            all_results.append(df)

        if dfs:
            results[str(cs)] = pd.concat(dfs, ignore_index=True)
            
    if all_results:
        final_DE_table = pd.concat(all_results, ignore_index=True)
    else:
        final_DE_table = pd.DataFrame()

    final_path = out_base / "final_DE_table.csv"
    final_DE_table.to_csv(final_path, index=False)

    print(f"[INFO] Wrote final table to {final_path}")
    return results

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("data_path", help="Path to .h5ad")
    p.add_argument("--out-dir", default="./wilcoxon_results")

    # columns in adata.obs
    p.add_argument("--cell-state-col", required=True)
    p.add_argument("--condition-col", required=True)
    p.add_argument("--sample-id-col", required=True)

    # reference / control
    p.add_argument("--reference-value", required=True, help='e.g. "TSteady"')

    # optional
    p.add_argument("--cell-states-to-exclude", nargs="*", default=[])

    # cell type meta
    p.add_argument("--cell-type-col", type=str, required=True)
    p.add_argument("--cell-type-val", type=str, required=True)

    # kept for symmetry with preprocessing (currently unused in wilcoxon.py)
    p.add_argument("--threshold", type=float, default=0.0125)
    p.add_argument("--min-cells", type = int, default = 3)
    p.add_argument("--region", required=True)
    p.add_argument("--annotation-level", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--paper", required=True)
    p.add_argument("--year", required=True)

    args = p.parse_args()


    cfg = WilcoxonConfig(
        data_path=args.data_path,
        out_dir=args.out_dir,
        cell_state_col=args.cell_state_col,
        condition_col=args.condition_col,
        sample_id_col=args.sample_id_col,
        cell_states_to_exclude=args.cell_states_to_exclude,
        cell_type_col=args.cell_type_col,
        cell_type_value=args.cell_type_val,
        threshold=args.threshold,
        min_cells = args.min_cells,        
        region=args.region,
        annotation_level=args.annotation_level,
        species=args.species,
        paper=args.paper,
        year=args.year,
    )

    results = run_Wilcoxon(cfg)