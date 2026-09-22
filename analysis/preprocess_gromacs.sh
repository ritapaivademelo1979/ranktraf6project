#!/usr/bin/env bash
# Make a whole-molecule, centred, protein-only trajectory + matching .tpr for analysis.
#
# Usage: ./preprocess_gromacs.sh md.tpr md.xtc out_prefix [dt_ps]   (dt_ps: keep one frame every dt_ps)
#   -> out_prefix_complex.tpr, out_prefix_complex.xtc
#
# The complex-only trajectory is ~10-50x smaller than the solvated one, which
# makes the Python analysis much faster. Default group is "Protein" (TRAF6 +
# peptide, both are protein chains). If the peptide carries non-standard
# residues/caps that GROMACS does not count as Protein, make an index group
# first (gmx make_ndx) and set GROUP and NDX below.
set -euo pipefail

TPR=$1; XTC=$2; OUT=$3; DT=${4:-0}
GMX=${GMX:-gmx}
GROUP=${GROUP:-Protein}
NDX=${NDX:-}
ndx_arg=(); [[ -n "$NDX" ]] && ndx_arg=(-n "$NDX")
dt_arg=();  [[ "$DT" != "0" ]] && dt_arg=(-dt "$DT")

# 1) Complex-only topology for MDAnalysis (keeps bonds and charges).
echo "$GROUP" | $GMX convert-tpr -s "$TPR" "${ndx_arg[@]}" -o "${OUT}_complex.tpr"

# 2) Make molecules whole, put all atoms in the box.
echo "$GROUP" | $GMX trjconv -s "$TPR" -f "$XTC" "${ndx_arg[@]}" "${dt_arg[@]}" \
    -pbc mol -o "${OUT}_tmp_whole.xtc"

# 3) Centre on the complex and keep molecules whole (compact unit cell).
#    Group order: centring group, then output group.
printf '%s\n%s\n' "$GROUP" "$GROUP" | $GMX trjconv -s "${OUT}_complex.tpr" -f "${OUT}_tmp_whole.xtc" \
    -center -pbc mol -ur compact -o "${OUT}_complex.xtc"

rm -f "${OUT}_tmp_whole.xtc"
echo "Wrote ${OUT}_complex.tpr and ${OUT}_complex.xtc"
# Fitting is done inside the Python analysis, and the peptide is re-imaged next
# to TRAF6 on every frame, so no -fit or -pbc nojump step is needed here.
