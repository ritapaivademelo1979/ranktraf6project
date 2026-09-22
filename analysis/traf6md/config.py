from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import yaml


@dataclass
class Replica:
    system: str
    name: str
    top: Path
    traj: Path


@dataclass
class Config:
    receptor_sel: str
    peptide_sel: str
    out_dir: Path
    replicas: list[Replica]
    stride: int = 1
    skip_ns: float = 0.0
    contact_cutoff_nm: float = 0.45
    saltbridge_cutoff_nm: float = 0.40
    hbonds: bool = True
    hotspots: dict[str, str] = field(default_factory=dict)
    per_system_sel: dict[str, dict[str, str]] = field(default_factory=dict)
    receptor_offset: int = 0
    peptide_offset: int = 0

    def selections_for(self, system: str) -> tuple[str, str]:
        over = self.per_system_sel.get(system, {})
        return over.get("receptor", self.receptor_sel), over.get("peptide", self.peptide_sel)


def load_config(path: str | Path) -> Config:
    path = Path(path)
    raw = yaml.safe_load(path.read_text())
    base = path.parent
    sel = raw.get("selections", {})
    analysis = raw.get("analysis", {})

    replicas, per_system_sel = [], {}
    for sname, sdef in raw["systems"].items():
        if "selections" in sdef:
            per_system_sel[sname] = sdef["selections"]
        for rname, rdef in sdef["replicas"].items():
            replicas.append(Replica(
                system=str(sname), name=str(rname),
                top=(base / rdef["top"]).resolve(),
                traj=(base / rdef["traj"]).resolve(),
            ))

    return Config(
        receptor_sel=sel["receptor"],
        peptide_sel=sel["peptide"],
        out_dir=(base / raw.get("out_dir", "results")).resolve(),
        replicas=replicas,
        stride=int(analysis.get("stride", 1)),
        skip_ns=float(analysis.get("skip_ns", 0.0)),
        contact_cutoff_nm=float(analysis.get("contact_cutoff_nm", 0.45)),
        saltbridge_cutoff_nm=float(analysis.get("saltbridge_cutoff_nm", 0.40)),
        hbonds=bool(analysis.get("hbonds", True)),
        hotspots=dict(raw.get("hotspots") or {}),
        per_system_sel=per_system_sel,
        receptor_offset=int((raw.get("numbering") or {}).get("receptor_offset", 0)),
        peptide_offset=int((raw.get("numbering") or {}).get("peptide_offset", 0)),
    )
