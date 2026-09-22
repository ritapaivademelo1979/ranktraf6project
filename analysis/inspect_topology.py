#!/usr/bin/env python3
"""Print the segments/chains and residue ranges in a topology to help write the selections in config.yaml."""
from __future__ import annotations
import argparse
import MDAnalysis as mda


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("top", help="Topology (.tpr recommended; .gro/.pdb also work)")
    p.add_argument("--test-sel", nargs="*", default=[], help="Selection strings to test (prints atom/residue counts)")
    args = p.parse_args()

    u = mda.Universe(args.top)
    print(f"{u.atoms.n_atoms} atoms, {u.residues.n_residues} residues\n")
    print(f"{'segid':<32} {'chainID':<8} {'residues':>8}  first → last")
    for seg in u.segments:
        res = seg.residues
        chains = sorted(set(getattr(seg.atoms, "chainIDs", ["-"])))
        print(f"{seg.segid:<32} {','.join(chains)[:8]:<8} {res.n_residues:>8}  "
              f"{res[0].resname}{res[0].resid} → {res[-1].resname}{res[-1].resid}")
    for sel in args.test_sel:
        ag = u.select_atoms(sel)
        print(f"\n'{sel}': {ag.n_atoms} atoms, {ag.residues.n_residues} residues")
        if ag.n_atoms:
            print("  " + " ".join(f"{r.resname}{r.resid}" for r in ag.residues[:40])
                  + (" ..." if ag.residues.n_residues > 40 else ""))


if __name__ == "__main__":
    main()
