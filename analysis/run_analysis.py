#!/usr/bin/env python3
"""Analyse every replica listed in config.yaml. Writes CSVs to <out_dir>/<system>/<replica>/."""
from __future__ import annotations
import argparse
import time

from traf6md.config import load_config
from traf6md.metrics import analyse_replica


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("config", help="YAML config (see config.example.yaml)")
    p.add_argument("--only", nargs="*", default=None, help="Analyse only these systems")
    p.add_argument("--force", action="store_true", help="Re-run replicas that already have results")
    args = p.parse_args()

    cfg = load_config(args.config)
    for rep in cfg.replicas:
        if args.only and rep.system not in args.only:
            continue
        out = cfg.out_dir / rep.system / rep.name
        if (out / "timeseries.csv").exists() and not args.force:
            print(f"[skip] {rep.system}/{rep.name} (already done; --force to redo)")
            continue
        t0 = time.time()
        print(f"[run ] {rep.system}/{rep.name}: {rep.traj.name}")
        results = analyse_replica(rep, cfg)
        out.mkdir(parents=True, exist_ok=True)
        for name, df in results.items():
            df.to_csv(out / f"{name}.csv", index=False, float_format="%.5g")
        print(f"[done] {rep.system}/{rep.name}: {len(results['timeseries'])} frames in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
