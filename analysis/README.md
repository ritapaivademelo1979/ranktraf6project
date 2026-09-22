# MD analysis: peptide–TRAF6 complexes

Scripts to analyse GROMACS trajectories of peptide–TRAF6 complexes and compare
or rank peptides across replicas.

| Output | Contents |
|---|---|
| `timeseries.csv` | TRAF6 backbone RMSD, peptide pose RMSD (after fitting on TRAF6), peptide conformational RMSD, Rg, COM and minimum distance, heavy-atom contacts, native-contact fraction Q, H-bonds, salt bridges |
| `rmsf_*.csv` | Cα RMSF for TRAF6 and the peptide |
| `contacts_*` | per-residue and residue-pair contact frequencies |
| `hbonds.csv`, `saltbridges.csv` | interface interaction occupancies |
| `summary/ranking.csv` | per-system mean ± SD over replicas, ranked by metric |
| `figures/` | time series, RMSF, contact heatmap, per-system metrics (PNG + SVG) |

## Setup

```bash
pip install -r analysis/requirements.txt
```

## The trajectory (.xtc) is too big for GitHub

You don't need the full solvated trajectory. Two options:

**A. Reduce it on the cluster, then copy or commit the small file (recommended).**
This keeps only the protein atoms (peptide + TRAF6, about 2,900 of the 38,000 atoms), makes molecules whole, and
optionally writes one frame every 100 ps:

```bash
cd analysis
./preprocess_gromacs.sh PEP2A_TRAF6_md.tpr PEP2A_md.xtc ../trajectories/PEP2A_r1 100
# -> ../trajectories/PEP2A_r1_complex.xtc  and  PEP2A_r1_complex.tpr
```

A 400 ns run at 100 ps per frame gives about 4,000 frames × 2,900 atoms, which is
roughly 40–60 MB. That's under GitHub's 100 MB per-file limit. For full resolution
(10 ps), use Git LFS or deposit the trajectories on Zenodo.

If you use the reduced trajectory, point `top:` in `config.yaml` at the matching
`*_complex.tpr` that the script writes, because the full `.tpr` has a different atom count.

**B. Run the analysis where the trajectory lives** (the cluster), and commit only
`results/`. The CSVs are small, and `summarize.py` can be re-run anywhere.

## Run

```bash
cd analysis
python inspect_topology.py ../PEP2A_TRAF6_md.tpr      # check segments / selections
python run_analysis.py config.yaml                     # per replica → ../results/<system>/<replica>/
python summarize.py config.yaml                        # → ../results/summary, ../results/figures
```

To add peptides or replicas, add entries under `systems:` in `config.yaml`.
Replicas that already have results are skipped; use `--force` to redo them.

## Notes

- Distance-based metrics use minimum-image distances. RMSD, RMSF and Rg need
  whole molecules, which `preprocess_gromacs.sh` ensures. The peptide is also
  re-imaged next to TRAF6 on every frame.
- The reference structure for RMSD and native contacts is the first frame of each replica.
- H-bonds: donor–acceptor ≤ 0.35 nm, D–H–A ≥ 150°. Salt bridges: acid O to base N ≤ 0.40 nm.
- The ranking averages per-metric ranks. It describes stability and interface
  persistence, not binding free energy. Use MM/PBSA (e.g. gmx_MMPBSA) or
  free-energy methods for energetics.
- `tests/make_synthetic.py` builds a small fake dataset to check the pipeline.
