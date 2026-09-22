#!/usr/bin/env python3
"""Combine replica results into per-system summaries, a ranking table, and figures.

Reads <out_dir>/<system>/<replica>/*.csv written by run_analysis.py and writes
<out_dir>/summary/ (CSV tables) and <out_dir>/figures/ (PNG + SVG).
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from traf6md.config import load_config
from traf6md import plotting

# metric column -> (label, higher_is_better) used for the ranking
RANK_METRICS = {
    "rmsd_peptide_pose_nm": ("Peptide pose RMSD (nm)", False),
    "q_native": ("Native contacts Q", True),
    "bound_fraction": ("Bound fraction", True),
    "contacts": ("Heavy-atom contacts", True),
    "hbonds": ("Interface H-bonds", True),
    "saltbridges": ("Salt bridges", True),
}


def load_tables(out_dir: Path, name: str) -> pd.DataFrame:
    frames = []
    for f in sorted(out_dir.glob(f"*/*/{name}.csv")):
        if f.parts[-3] in ("summary", "figures"):
            continue
        df = pd.read_csv(f)
        df.insert(0, "replica", f.parent.name)
        df.insert(0, "system", f.parent.parent.name)
        frames.append(df)
    if not frames:
        raise SystemExit(f"No {name}.csv found under {out_dir}. Run run_analysis.py first.")
    return pd.concat(frames, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("config")
    p.add_argument("--smooth-ns", type=float, default=1.0, help="Running-mean window for time-series plots (ns)")
    p.add_argument("--heatmap-min", type=float, default=0.3,
                   help="Show receptor residues contacted at least this often in some system (0-1)")
    args = p.parse_args()

    cfg = load_config(args.config)
    out = cfg.out_dir
    sum_dir, fig_dir = out / "summary", out / "figures"
    sum_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    order = list(dict.fromkeys(r.system for r in cfg.replicas))  # config order = colour order

    ts = load_tables(out, "timeseries")
    order = [s for s in order if s in set(ts.system)]
    eq = ts[ts.time_ns >= cfg.skip_ns].copy()
    eq["bound_fraction"] = (eq.mindist_nm < cfg.contact_cutoff_nm).astype(float)

    metric_cols = [c for c in eq.columns if c not in ("system", "replica", "time_ns")]
    per_rep = eq.groupby(["system", "replica"])[metric_cols].mean().reset_index()
    per_rep["frames"] = eq.groupby(["system", "replica"]).size().values
    per_rep.to_csv(sum_dir / "per_replica.csv", index=False, float_format="%.4g")

    agg = per_rep.groupby("system")[metric_cols].agg(["mean", "std"])
    agg.columns = [f"{m}_{s}" for m, s in agg.columns]
    agg.insert(0, "n_replicas", per_rep.groupby("system").size())
    agg = agg.reindex(order)
    agg.to_csv(sum_dir / "per_system.csv", float_format="%.4g")

    # Ranking: rank each metric across systems (1 = best), then average the ranks.
    rank = pd.DataFrame(index=agg.index)
    for m, (_label, higher) in RANK_METRICS.items():
        if f"{m}_mean" in agg:
            rank[m] = agg[f"{m}_mean"]
            rank[f"rank_{m}"] = agg[f"{m}_mean"].rank(ascending=not higher, method="min")
    rank_cols = [c for c in rank.columns if c.startswith("rank_")]
    rank["mean_rank"] = rank[rank_cols].mean(axis=1)
    rank = rank.sort_values("mean_rank")
    rank.insert(0, "position", np.arange(1, len(rank) + 1))
    rank.to_csv(sum_dir / "ranking.csv", float_format="%.4g")

    # Residue-level tables averaged over replicas
    tables = {}
    for name, key in [("contacts_receptor_residues", "contact_frequency"),
                      ("contacts_peptide_residues", "contact_frequency"),
                      ("rmsf_receptor", "rmsf_nm"), ("rmsf_peptide", "rmsf_nm")]:
        df = load_tables(out, name)
        g = df.groupby(["system", "resid", "resname"])[key].agg(["mean", "std"]).reset_index()
        g.to_csv(sum_dir / f"{name}.csv", index=False, float_format="%.4g")
        tables[name] = g
    for name, key, cols in [("hbonds", "occupancy", ["donor_res", "donor_atom", "acceptor_res", "acceptor_atom"]),
                            ("saltbridges", "occupancy", ["acid_res", "base_res"]),
                            ("contacts_residue_pairs", "frequency", ["receptor_res", "peptide_res"])]:
        try:
            df = load_tables(out, name)
        except SystemExit:
            continue
        n_rep = per_rep.groupby("system").size()
        # Missing interactions in a replica count as 0 occupancy, so divide by the replica count.
        g = df.groupby(["system", *cols])[key].sum().div(n_rep, level="system").rename(f"mean_{key}")
        g.reset_index().sort_values(["system", f"mean_{key}"], ascending=[True, False]) \
            .to_csv(sum_dir / f"{name}.csv", index=False, float_format="%.4g")

    plotting.timeseries(ts, order, fig_dir, cfg.skip_ns, args.smooth_ns)
    plotting.metric_dots(per_rep, order, RANK_METRICS, fig_dir)
    plotting.rmsf(tables["rmsf_receptor"], order, fig_dir / "rmsf_receptor", "TRAF6 residue")
    plotting.rmsf(tables["rmsf_peptide"], order, fig_dir / "rmsf_peptide", "Peptide residue")
    plotting.contact_heatmap(tables["contacts_receptor_residues"], order, fig_dir, args.heatmap_min)

    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print("\nRanking (equilibrated window t >= %.0f ns; 1 = best):" % cfg.skip_ns)
        print(rank[["position", "mean_rank", *[m for m in RANK_METRICS if m in rank]]].round(3).to_string())
    print(f"\nTables  → {sum_dir}\nFigures → {fig_dir}")


if __name__ == "__main__":
    main()
