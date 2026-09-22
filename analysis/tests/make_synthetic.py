#!/usr/bin/env python3
"""Build a small fake receptor–peptide trajectory (PDB + XTC) to smoke-test the pipeline.

The peptide drifts away in replica 'r2' of system 'weak' and is moved by one box
vector in some frames to check that periodic imaging is handled.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import MDAnalysis as mda

RES_ATOMS = ["N", "H", "CA", "C", "O"]
SIDE = {"ARG": ["CB", "CZ", "NH1", "NH2"], "GLU": ["CB", "CD", "OE1", "OE2"],
        "LYS": ["CB", "NZ"], "ALA": ["CB"], "PRO": ["CB"]}


def build(rec_seq, pep_seq):
    names, resn, resid_of_atom, seg_of_res, resids, resnames = [], [], [], [], [], []
    for seg, seq in enumerate([rec_seq, pep_seq]):
        for i, rn in enumerate(seq):
            r = len(resids)
            resids.append(i + 1); resnames.append(rn); seg_of_res.append(seg)
            for a in RES_ATOMS + SIDE[rn]:
                names.append(a); resid_of_atom.append(r)
    u = mda.Universe.empty(len(names), n_residues=len(resids), n_segments=2,
                           atom_resindex=resid_of_atom, residue_segindex=seg_of_res, trajectory=True)
    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("type", [n[0] for n in names])
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("resid", resids)
    u.add_TopologyAttr("segid", ["REC", "PEP"])
    u.add_TopologyAttr("chainIDs", [ "A" if seg_of_res[r] == 0 else "B" for r in resid_of_atom])
    return u


def coords(u, rng):
    x = np.zeros((u.atoms.n_atoms, 3))
    for res in u.residues:
        base = (np.array([3.8 * res.resid, 0, 0]) if res.segment.segid == "REC"
                else np.array([3.8 * res.resid + 20, 4.5, 0]))
        for k, a in enumerate(res.atoms):
            x[a.index] = base + rng.normal(0, 0.8, 3) + np.array([0, 0.6 * k * (1 if res.segment.segid == "REC" else 0.3), 0])
    # put the peptide ARG/GLU side-chain tips next to each other to form a salt bridge
    return x


def main(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    rec_seq = ["ALA", "ARG", "LYS", "GLU", "ALA", "ARG", "PRO", "ALA", "GLU", "LYS", "ALA", "ARG"] * 2
    pep_seq = ["PRO", "GLU", "GLU", "ALA", "ARG", "GLU"]
    u = build(rec_seq, pep_seq)
    base = coords(u, rng)
    pep = u.select_atoms("segid PEP")
    box = np.array([120.0, 120.0, 120.0, 90, 90, 90])
    for system, drift in [("strong", [0.0, 0.0]), ("weak", [0.0, 0.25])]:
        for r, d in enumerate(drift, start=1):
            top = out / f"{system}_r{r}.pdb"
            traj = out / f"{system}_r{r}.xtc"
            u.atoms.positions = base
            u.dimensions = box
            u.atoms.write(top)
            with mda.Writer(str(traj), u.atoms.n_atoms) as w:
                for f in range(200):
                    x = base + rng.normal(0, 0.4, base.shape)
                    x[pep.indices] += np.array([0, 0, d * f])
                    if f % 7 == 3:  # jump the peptide to a neighbouring periodic image
                        x[pep.indices] += np.array([box[0], 0, 0])
                    u.atoms.positions = x
                    u.dimensions = box
                    u.trajectory.ts.time = f * 100.0  # 0.1 ns per frame
                    w.write(u.atoms)
    cfg = out / "config.yaml"
    cfg.write_text(f"""out_dir: results
selections: {{receptor: "segid REC", peptide: "segid PEP"}}
analysis: {{stride: 1, skip_ns: 2, contact_cutoff_nm: 0.45, hbonds: true}}
hotspots: {{basic_patch: "resname ARG LYS and resid 1-6"}}
systems:
  strong:
    replicas:
      r1: {{top: strong_r1.pdb, traj: strong_r1.xtc}}
      r2: {{top: strong_r2.pdb, traj: strong_r2.xtc}}
  weak:
    replicas:
      r1: {{top: weak_r1.pdb, traj: weak_r1.xtc}}
      r2: {{top: weak_r2.pdb, traj: weak_r2.xtc}}
""")
    print(f"Synthetic data and config written to {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "synthetic"))
