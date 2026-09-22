"""Per-replica analysis of a peptide–TRAF6 trajectory.

All distance-based quantities (contacts, minimum distances, salt bridges,
H-bonds) use minimum-image distances, so they don't depend on how the
trajectory was wrapped. RMSD, RMSF and pose-dependent quantities also need the
receptor and peptide to be whole molecules. `preprocess_gromacs.sh` takes care
of that. The peptide is moved into the periodic image nearest the receptor on
every frame before fitting.

MDAnalysis works in Angstrom, and every output is converted to nm.
"""
from __future__ import annotations

from collections import Counter
import warnings

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import distance_array, minimize_vectors

from .config import Config, Replica

# AMBER and CHARMM names (C-terminal oxygens: OC1/OC2 in AMBER, OT1/OT2 in CHARMM; protonated His: HIP/HSP).
ACIDIC = ("((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)"
          " or name OC1 OC2 OT1 OT2 OXT)")
BASIC = ("((resname ARG and name NH1 NH2 NE) or (resname LYS and name NZ)"
         " or (resname HIP HSP and name ND1 NE2))")


# Residue-number offsets for output tables (set from the config in analyse_replica).
_OFFSET = {"receptor": 0, "peptide": 0}


def _res_label(res, peptide: bool = False) -> str:
    """Residue label; peptide residues get a 'pep:' prefix so numbering can't clash with the receptor."""
    off = _OFFSET["peptide" if peptide else "receptor"]
    return f"{'pep:' if peptide else ''}{res.resname}{res.resid + off}"


def _heavy(ag: mda.AtomGroup) -> mda.AtomGroup:
    return ag.select_atoms("not name H*")


def _check_whole(ag: mda.AtomGroup, label: str) -> None:
    """Warn if bonded atoms are > 0.4 nm apart (molecule split across the box)."""
    try:
        bonds = ag.intra_bonds
    except Exception:  # no bond information (e.g. PDB topology)
        return
    if len(bonds) == 0:
        return
    if bonds.values().max() > 4.0:
        warnings.warn(f"{label} is not whole in frame 0 (max bond {bonds.values().max() / 10:.2f} nm). "
                      "Run preprocess_gromacs.sh first, or RMSD/RMSF/Rg will be wrong.")


class _ResidueIndex:
    """Map atoms of a group to 0..n_residues-1 of that group."""

    def __init__(self, group: mda.AtomGroup):
        self.residues = group.residues
        lut = {ri: k for k, ri in enumerate(self.residues.resindices)}
        self.lut = lut

    def of(self, atoms: mda.AtomGroup) -> np.ndarray:
        return np.array([self.lut[ri] for ri in atoms.resindices], dtype=int)


def analyse_replica(rep: Replica, cfg: Config) -> dict[str, pd.DataFrame]:
    rec_sel, pep_sel = cfg.selections_for(rep.system)
    _OFFSET["receptor"], _OFFSET["peptide"] = cfg.receptor_offset, cfg.peptide_offset
    u = mda.Universe(str(rep.top), str(rep.traj))
    rec = u.select_atoms(rec_sel)
    pep = u.select_atoms(pep_sel)
    if rec.n_atoms == 0 or pep.n_atoms == 0:
        raise ValueError(f"[{rep.system}/{rep.name}] empty selection: receptor={rec.n_atoms} "
                         f"peptide={pep.n_atoms} atoms. Use inspect_topology.py to find the right selections.")
    if len(rec & pep):
        raise ValueError(f"[{rep.system}/{rep.name}] receptor and peptide selections overlap.")

    rec_bb, pep_bb = rec.select_atoms("backbone"), pep.select_atoms("backbone")
    rec_ca, pep_ca = rec.select_atoms("name CA"), pep.select_atoms("name CA")
    rec_h, pep_h = _heavy(rec), _heavy(pep)
    rec_idx, pep_idx = _ResidueIndex(rec), _ResidueIndex(pep)
    rec_h_res, pep_h_res = rec_idx.of(rec_h), pep_idx.of(pep_h)

    hotspots = {}
    for name, sel in cfg.hotspots.items():
        ag = _heavy(rec.select_atoms(sel))
        if ag.n_atoms == 0:
            warnings.warn(f"hotspot '{name}' ({sel}) matched no receptor atoms, skipping")
            continue
        hotspots[name] = ag

    # Salt bridges in both directions: peptide acid–receptor base and peptide base–receptor acid.
    sb_sets = [(pep.select_atoms(ACIDIC), rec.select_atoms(BASIC)),
               (pep.select_atoms(BASIC), rec.select_atoms(ACIDIC))]
    sb_sets = [(a, b) for a, b in sb_sets if a.n_atoms and b.n_atoms]

    pep_ix = set(pep.indices)
    cut_A = cfg.contact_cutoff_nm * 10.0
    sb_cut_A = cfg.saltbridge_cutoff_nm * 10.0
    skip_ps = cfg.skip_ns * 1000.0

    # Reference is frame 0, with the peptide moved into the receptor's image.
    u.trajectory[0]
    _check_whole(rec, "receptor")
    _check_whole(pep, "peptide")
    shift0 = _image_shift(rec, pep, u.trajectory.ts.dimensions)
    ref_rec_bb = rec_bb.positions.copy()
    ref_center = ref_rec_bb.mean(axis=0)
    ref_rec_bb_c = ref_rec_bb - ref_center
    ref_pep_bb = pep_bb.positions + shift0
    d0 = distance_array(rec_h.positions, pep_h.positions, box=u.trajectory.ts.dimensions)
    native = set(zip(rec_h_res[np.nonzero(d0 < cut_A)[0]], pep_h_res[np.nonzero(d0 < cut_A)[1]]))
    if not native:
        warnings.warn(f"[{rep.system}/{rep.name}] no receptor–peptide contacts in frame 0; Q will be NaN.")

    rows = []
    n_eq = 0
    rec_ca_sum = np.zeros((rec_ca.n_atoms, 3)); rec_ca_sq = np.zeros(rec_ca.n_atoms)
    pep_ca_sum = np.zeros((pep_ca.n_atoms, 3)); pep_ca_sq = np.zeros(pep_ca.n_atoms)
    rec_res_hits = np.zeros(len(rec_idx.residues))
    pep_res_hits = np.zeros(len(pep_idx.residues))
    pair_hits: Counter = Counter()
    sb_hits: Counter = Counter()

    for ts in u.trajectory[::cfg.stride]:
        box = ts.dimensions
        shift = _image_shift(rec, pep, box)

        # Superpose on receptor backbone (applied to copies, raw coordinates stay put).
        mob = rec_bb.positions
        mob_center = mob.mean(axis=0)
        R, rmsd_rec = align.rotation_matrix(mob - mob_center, ref_rec_bb_c)

        def fit(x):
            return (x - mob_center) @ R.T + ref_center

        pep_bb_fit = fit(pep_bb.positions + shift)
        pose_rmsd = np.sqrt(((pep_bb_fit - ref_pep_bb) ** 2).sum(axis=1).mean())
        conf_rmsd = rms.rmsd(pep_bb.positions, ref_pep_bb, center=True, superposition=True)

        d = distance_array(rec_h.positions, pep_h.positions, box=box)
        ci, cj = np.nonzero(d < cut_A)
        res_pairs = set(zip(rec_h_res[ci], pep_h_res[cj]))
        com_vec = minimize_vectors(pep.center_of_mass() - rec.center_of_mass(), box)

        row = {
            "time_ns": ts.time / 1000.0,
            "rmsd_receptor_bb_nm": rmsd_rec / 10.0,
            "rmsd_peptide_pose_nm": pose_rmsd / 10.0,
            "rmsd_peptide_conf_nm": conf_rmsd / 10.0,
            "rg_receptor_nm": rec.radius_of_gyration() / 10.0,
            "rg_peptide_nm": pep.radius_of_gyration() / 10.0,
            "com_distance_nm": np.linalg.norm(com_vec) / 10.0,
            "mindist_nm": d.min() / 10.0,
            "contacts": len(ci),
            "residue_contacts": len(res_pairs),
            "q_native": len(res_pairs & native) / len(native) if native else np.nan,
        }
        for name, ag in hotspots.items():
            row[f"mindist_{name}_nm"] = distance_array(ag.positions, pep_h.positions, box=box).min() / 10.0

        frame_sb = set()
        for acid, base in sb_sets:
            dd = distance_array(acid.positions, base.positions, box=box)
            ai, bi = np.nonzero(dd < sb_cut_A)
            for a, b in zip(ai, bi):
                acid_on_pep = acid[a].index in pep_ix
                frame_sb.add((_res_label(acid[a].residue, acid_on_pep),
                              _res_label(base[b].residue, not acid_on_pep),
                              "peptide" if acid_on_pep else "receptor"))
        row["saltbridges"] = len(frame_sb)
        rows.append(row)

        if ts.time >= skip_ps:
            n_eq += 1
            x = fit(rec_ca.positions)
            rec_ca_sum += x; rec_ca_sq += (x ** 2).sum(axis=1)
            x = fit(pep_ca.positions + shift)
            pep_ca_sum += x; pep_ca_sq += (x ** 2).sum(axis=1)
            if len(ci):
                rec_res_hits[np.unique(rec_h_res[ci])] += 1
                pep_res_hits[np.unique(pep_h_res[cj])] += 1
            pair_hits.update(res_pairs)
            sb_hits.update(frame_sb)

    if n_eq == 0:
        raise ValueError(f"[{rep.system}/{rep.name}] no frames after skip_ns={cfg.skip_ns}.")

    out = {"timeseries": pd.DataFrame(rows)}
    ro, po = cfg.receptor_offset, cfg.peptide_offset
    out["rmsf_receptor"] = _rmsf_df(rec_ca, rec_ca_sum, rec_ca_sq, n_eq, ro)
    out["rmsf_peptide"] = _rmsf_df(pep_ca, pep_ca_sum, pep_ca_sq, n_eq, po)
    out["contacts_receptor_residues"] = _res_freq_df(rec_idx.residues, rec_res_hits / n_eq, ro)
    out["contacts_peptide_residues"] = _res_freq_df(pep_idx.residues, pep_res_hits / n_eq, po)
    out["contacts_residue_pairs"] = pd.DataFrame(
        [{"receptor_res": _res_label(rec_idx.residues[i]), "peptide_res": _res_label(pep_idx.residues[j], True),
          "native": (i, j) in native, "frequency": c / n_eq} for (i, j), c in pair_hits.items()],
        columns=["receptor_res", "peptide_res", "native", "frequency"],
    ).sort_values("frequency", ascending=False)
    out["saltbridges"] = pd.DataFrame(
        [{"acid_res": a, "base_res": b, "acid_on": side, "occupancy": c / n_eq}
         for (a, b, side), c in sb_hits.items()],
        columns=["acid_res", "base_res", "acid_on", "occupancy"],
    ).sort_values("occupancy", ascending=False)

    if cfg.hbonds:
        hb_ts, hb_occ = _hbonds(u, rec_sel, pep_sel, cfg.stride, skip_ps)
        out["timeseries"]["hbonds"] = hb_ts.reindex(range(len(rows)), fill_value=0).to_numpy()
        out["hbonds"] = hb_occ
    return out


def _image_shift(rec, pep, box) -> np.ndarray:
    """Translation that puts the (whole) peptide in the periodic image closest to the receptor."""
    if box is None or not np.all(box[:3] > 0):
        return np.zeros(3)
    d = pep.center_of_geometry() - rec.center_of_geometry()
    return minimize_vectors(d, box) - d


def _rmsf_df(ca, s, sq, n, offset=0) -> pd.DataFrame:
    mean = s / n
    rmsf = np.sqrt(np.maximum(sq / n - (mean ** 2).sum(axis=1), 0.0))
    return pd.DataFrame({"resid": ca.resids + offset, "resname": ca.resnames, "rmsf_nm": rmsf / 10.0})


def _res_freq_df(residues, freq, offset=0) -> pd.DataFrame:
    return pd.DataFrame({"resid": residues.resids + offset, "resname": residues.resnames, "contact_frequency": freq})


def _hbonds(u, rec_sel, pep_sel, stride, skip_ps):
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

    both = f"({rec_sel}) or ({pep_sel})"
    if hasattr(u.atoms, "charges"):
        guesser = HydrogenBondAnalysis(u)
        h_sel, a_sel = guesser.guess_hydrogens(both), guesser.guess_acceptors(both)
    else:  # no charges (e.g. PDB topology): fall back to atom names
        h_sel, a_sel = f"({both}) and name H*", f"({both}) and (name O* or name N*)"
    # Without bonds, donor–H pairs are assigned by distance and need an explicit donor selection.
    d_sel = None if hasattr(u, "bonds") and len(u.bonds) else f"({both}) and (name O* or name N*)"
    hba = HydrogenBondAnalysis(u, between=[rec_sel, pep_sel], donors_sel=d_sel, hydrogens_sel=h_sel,
                               acceptors_sel=a_sel, d_a_cutoff=3.5, d_h_a_angle_cutoff=150)
    hba.run(step=stride)

    per_frame = pd.Series(hba.count_by_time(), dtype=int)
    hb = hba.results.hbonds
    times = dict(zip(hba.frames, hba.times))
    eq_frames = [f for f in hba.frames if times[f] >= skip_ps]
    rec_ix = set(u.select_atoms(rec_sel).indices)
    occ = Counter()
    for frame, donor, _h, acc, *_ in hb:
        if times[int(frame)] < skip_ps:
            continue
        D, A = u.atoms[int(donor)], u.atoms[int(acc)]
        d_rec, a_rec = D.index in rec_ix, A.index in rec_ix
        occ[(_res_label(D.residue, not d_rec), D.name, _res_label(A.residue, not a_rec), A.name,
             "receptor" if d_rec else "peptide")] += 1
    df = pd.DataFrame(
        [{"donor_res": dr, "donor_atom": da, "acceptor_res": ar, "acceptor_atom": aa,
          "donor_on": side, "occupancy": c / max(len(eq_frames), 1)}
         for (dr, da, ar, aa, side), c in occ.items()],
        columns=["donor_res", "donor_atom", "acceptor_res", "acceptor_atom", "donor_on", "occupancy"],
    ).sort_values("occupancy", ascending=False)
    return per_frame.reset_index(drop=True), df
