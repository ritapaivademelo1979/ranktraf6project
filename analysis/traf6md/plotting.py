"""Publication figures. One measure per axis, categorical colours in a fixed order per system."""
from __future__ import annotations
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# Validated colorblind-safe categorical order (adjacent-pair CVD ΔE ≥ 8 on a light background).
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
TEXT, TEXT_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("blue", ["#f7fbff", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])

TS_PANELS = [
    ("rmsd_receptor_bb_nm", "TRAF6 bb RMSD (nm)"),
    ("rmsd_peptide_pose_nm", "Pep. pose RMSD (nm)"),
    ("rmsd_peptide_conf_nm", "Pep. conf. RMSD (nm)"),
    ("q_native", "Native contacts Q"),
    ("contacts", "Heavy-atom contacts"),
    ("hbonds", "Interface H-bonds"),
    ("saltbridges", "Salt bridges"),
    ("com_distance_nm", "COM distance (nm)"),
    ("mindist_nm", "Min. distance (nm)"),
    ("rg_peptide_nm", "Peptide Rg (nm)"),
]

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT, "xtick.color": TEXT_2,
    "ytick.color": TEXT_2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
    "figure.dpi": 150, "savefig.bbox": "tight",
})


def colors_for(order: list[str]) -> dict[str, str]:
    if len(order) > len(PALETTE):
        warnings.warn(f"{len(order)} systems but only {len(PALETTE)} distinguishable colours; "
                      "colours repeat, so split the systems across several configs or read the summary tables.")
    return {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(order)}


def _save(fig, path: Path):
    for ext in ("png", "svg"):
        fig.savefig(path.with_suffix(f".{ext}"))
    plt.close(fig)


def timeseries(ts: pd.DataFrame, order, fig_dir: Path, skip_ns: float, smooth_ns: float):
    col = colors_for(order)
    panels = [(c, l) for c, l in TS_PANELS if c in ts] + \
             [(c, c.replace("mindist_", "Min. dist. ").replace("_nm", " (nm)"))
              for c in ts.columns if c.startswith("mindist_") and c != "mindist_nm"]
    ncol = 2
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 2.0 * nrow), sharex=True, squeeze=False)
    for ax, (c, label) in zip(axes.flat, panels):
        for s in order:
            d = ts[ts.system == s]
            wide = d.pivot_table(index="time_ns", columns="replica", values=c)
            dt = np.median(np.diff(wide.index)) if len(wide) > 1 else 1.0
            win = max(int(round(smooth_ns / dt)), 1)
            wide = wide.rolling(win, center=True, min_periods=1).mean()
            m, sd = wide.mean(axis=1), wide.std(axis=1)
            ax.plot(m.index, m, color=col[s], lw=1.2, label=s)
            if wide.shape[1] > 1:
                ax.fill_between(m.index, m - sd, m + sd, color=col[s], alpha=0.15, lw=0)
        if skip_ns > 0:
            ax.axvspan(ts.time_ns.min(), skip_ns, color=GRID, alpha=0.5, lw=0, zorder=0)
        ax.set_ylabel(label)
    for ax in axes.flat[len(panels):]:
        ax.remove()
    for ax in axes.flat[max(len(panels) - ncol, 0):len(panels)]:
        ax.set_xlabel("Time (ns)")
        ax.xaxis.set_tick_params(labelbottom=True)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=min(len(order), 4), bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Replica mean ± SD (shaded); grey = equilibration window excluded from statistics",
                 y=1.05, fontsize=8, color=TEXT_2)
    live = list(axes.flat[:len(panels)])
    live[0].set_xlim(ts.time_ns.min(), ts.time_ns.max())
    fig.tight_layout()
    fig.align_ylabels(live[0::2])
    fig.align_ylabels(live[1::2])
    _save(fig, fig_dir / "timeseries")


def metric_dots(per_rep: pd.DataFrame, order, metrics: dict, fig_dir: Path):
    """Per-system mean ± SD across replicas with individual replicas as points."""
    cols = [m for m in metrics if m in per_rep]
    fig, axes = plt.subplots(1, len(cols), figsize=(2.0 * len(cols), 0.35 * len(order) + 1.0), sharey=True)
    axes = np.atleast_1d(axes)
    y = np.arange(len(order))[::-1]
    for ax, m in zip(axes, cols):
        for yi, s in zip(y, order):
            v = per_rep.loc[per_rep.system == s, m].to_numpy()
            ax.scatter(v, np.full(len(v), yi), s=14, color=PALETTE[0], alpha=0.45, lw=0)
            ax.errorbar(v.mean(), yi, xerr=v.std(ddof=1) if len(v) > 1 else None, fmt="o", ms=5,
                        color=PALETTE[0], ecolor=TEXT_2, elinewidth=1, capsize=2)
        label, higher = metrics[m]
        ax.set_title(f"{label}\n({'higher' if higher else 'lower'} = better)", fontsize=8, color=TEXT)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(y, order)
    fig.tight_layout()
    _save(fig, fig_dir / "metrics_by_system")


def rmsf(df: pd.DataFrame, order, path: Path, xlabel: str):
    col = colors_for(order)
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    for s in order:
        d = df[df.system == s].sort_values("resid")
        ax.plot(d.resid, d["mean"], color=col[s], lw=1.2, label=s)
        ax.fill_between(d.resid, d["mean"] - d["std"].fillna(0), d["mean"] + d["std"].fillna(0),
                        color=col[s], alpha=0.15, lw=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cα RMSF (nm)")
    if len(order) > 1:
        ax.legend(ncol=min(len(order), 4), fontsize=8)
    fig.tight_layout()
    _save(fig, path)


def contact_heatmap(df: pd.DataFrame, order, fig_dir: Path, min_freq: float):
    wide = df.pivot_table(index="system", columns=["resid", "resname"], values="mean").reindex(order)
    keep = wide.columns[(wide.max(axis=0) >= min_freq).to_numpy()]
    if len(keep) == 0:
        warnings.warn(f"No receptor residue reaches contact frequency {min_freq}; heatmap skipped.")
        return
    wide = wide[keep]
    fig, ax = plt.subplots(figsize=(max(3.0, 0.22 * len(keep) + 1.5), 0.35 * len(order) + 1.2))
    im = ax.imshow(wide.to_numpy(), aspect="auto", cmap=SEQ, vmin=0, vmax=1)
    ax.set_xticks(range(len(keep)), [f"{rn}{ri}" for ri, rn in keep], rotation=90, fontsize=7)
    ax.set_yticks(range(len(order)), order)
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Contact frequency")
    ax.set_title(f"TRAF6 residues contacting the peptide\n(≥ {min_freq:.0%} of frames in any system)",
                 fontsize=8, color=TEXT)
    fig.tight_layout()
    _save(fig, fig_dir / "contact_heatmap_receptor")
