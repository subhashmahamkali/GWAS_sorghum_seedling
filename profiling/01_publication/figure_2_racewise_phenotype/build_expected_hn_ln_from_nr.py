#!/usr/bin/env python3
"""Build HN/LN expected-value BLUE table from existing BLUE + NR columns.

This script does not modify the source file. It writes:
1) A new BLUE table where HN/LN columns are replaced with expected values
   (BLUE + intercept recovered per year/trait/treatment).
2) A small intercept summary table for auditability.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path("/Users/subhashmahamkali/Documents/gwas_sap")
SRC = REPO_ROOT / "data/1.Phenotype_data/1.2020_2021_SAP/1.BLUEs_SAP_2020_2021.csv"
OUT = REPO_ROOT / "data/1.Phenotype_data/1.2020_2021_SAP/1.BLUEs_SAP_2020_2021_HN_LN_expected.csv"
OUT_INTERCEPTS = (
    REPO_ROOT
    / "data/1.Phenotype_data/1.2020_2021_SAP/1.BLUEs_SAP_2020_2021_HN_LN_expected_intercepts.csv"
)


def main() -> None:
    df = pd.read_csv(SRC)
    out_df = df.copy()

    nr_pattern = re.compile(r"^(\d{4}):(.+):NR$")
    nr_cols = [c for c in df.columns if nr_pattern.match(c)]
    rows = []

    for nr_col in nr_cols:
        year, trait = nr_pattern.match(nr_col).groups()
        hn_col = f"{year}:HN:{trait}"
        ln_col = f"{year}:LN:{trait}"

        if hn_col not in df.columns or ln_col not in df.columns:
            continue

        fit = df[[hn_col, ln_col, nr_col]].dropna()
        if len(fit) < 3:
            continue

        bh = fit[hn_col].to_numpy(dtype=float)
        bl = fit[ln_col].to_numpy(dtype=float)
        nr = fit[nr_col].to_numpy(dtype=float)

        # Rearranged from:
        # nr = ((bh + i_hn) - (bl + i_ln)) / (bl + i_ln)
        # -> i_hn - (1 + nr) * i_ln = (1 + nr) * bl - bh
        a = np.column_stack([np.ones(len(fit)), -(1.0 + nr)])
        y = (1.0 + nr) * bl - bh
        intercept_hn, intercept_ln = np.linalg.lstsq(a, y, rcond=None)[0]

        out_df[hn_col] = df[hn_col] + intercept_hn
        out_df[ln_col] = df[ln_col] + intercept_ln

        pred_nr = ((bh + intercept_hn) - (bl + intercept_ln)) / (bl + intercept_ln)
        rmse = float(np.sqrt(np.mean((pred_nr - nr) ** 2)))

        rows.append(
            {
                "year": year,
                "trait": trait,
                "n_fitted": int(len(fit)),
                "intercept_hn": float(intercept_hn),
                "intercept_ln": float(intercept_ln),
                "nr_reconstruction_rmse": rmse,
            }
        )

    intercepts = pd.DataFrame(rows).sort_values(["year", "trait"]).reset_index(drop=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT, index=False)
    intercepts.to_csv(OUT_INTERCEPTS, index=False)

    print(f"Wrote expected BLUE table: {OUT}")
    print(f"Wrote intercept summary: {OUT_INTERCEPTS}")
    if not intercepts.empty:
        print(
            "Intercept recovery RMSE (median/max): "
            f"{intercepts['nr_reconstruction_rmse'].median():.3e} / "
            f"{intercepts['nr_reconstruction_rmse'].max():.3e}"
        )


if __name__ == "__main__":
    main()
