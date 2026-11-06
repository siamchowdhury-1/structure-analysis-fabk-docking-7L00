<!--
============================================================
Siam Chowdhury (Beck Scholar); [www.siamchowdhury.com]
Computational and Medicinal Chemistry
Mentor: Mohammad Alam PhD., Professor of Chemistry
[Dr. Alam's Research Team] www.alamresearch.org
Beck College of Sciences and Mathematics, Arkansas State University, AR, USA
----------------------------------------------------------------------------
-->

# Structure Analysis of the Interaction of 7L00 (FabK) & Compounds (5sa23, 17sa23, 20sa23, 22sa23, 25sa23)

PyMOL scripts for analysing and rendering the binding interactions between the
enoyl-ACP reductase FabK (PDB 7L00) and the s_series compounds.

**Siam Chowdhury** (Beck Scholar) · [www.siamchowdhury.com](https://www.siamchowdhury.com)
Computational and Medicinal Chemistry
Mentor: **Mohammad Alam, PhD**, Professor of Chemistry
[Dr. Alam's Research Team](https://www.alamresearch.org)
Beck College of Sciences and Mathematics, Arkansas State University, AR, USA

---

## Contents

| File | Purpose |
|---|---|
| `make_docking_figures.py` | Binding-pose overview, close-up with H-bonds, and an all-compound overlay |
| `make_interaction_figures.py` | Interaction map: H-bonds, halogen bonds, salt bridges, π–π, π–cation, hydrophobic |
| `acs_labels.pml` | ACS-style labelling of a single compound (NH, OH, NH₂ …) |

---

## Requirements

PyMOL (open-source or Incentive) with Python 3. Everything else ships with PyMOL
(`chempy.cpv` is used for the ring geometry).

```bash
conda activate docking
pymol -cq --version
```

---

## Input files expected

> **Note on paths.** These scripts were developed and run in:
>
> ```
> ~/Documents/Computational_work/WORK LAB/s_series
> ```
>
> Nothing is hard-coded to that location — every path in the scripts is relative
> to the current working directory. Place the receptor, crystal ligand and
> `docking_*` folders together as shown below, `cd` into that folder, and run
> from there. If your layout differs, edit the `RECEPTOR`, `CRYSTAL` and
> `COMPOUNDS` constants at the top of each script.

```
s_series/
├── 7L00_AB_protein_FMN.pdb        # receptor: protein + FMN cofactor
├── XCJ_A402_crystal.sdf           # crystal ligand (sets the common camera view)
├── docking_focused_<cmpd>/
│   └── <cmpd>_best.mol2           # top-scoring pose, read by both scripts
└── docking_<cmpd>/                # same layout, non-focused run
```

`*_best.mol2` is used rather than `.pdbqt` because it preserves bond orders,
which the aromatic-ring and valence detection depend on.

---

# Workflow

This covers the visualisation and interaction-analysis stage only. Receptor
preparation, redocking validation, Vina docking and MM-GBSA are upstream steps
documented separately.

Prerequisite: docking already complete, with a top-scoring pose per compound at
`docking_focused_<cmpd>/<cmpd>_best.mol2`.

```bash
conda activate docking
cd ~/Documents/Computational_work/"WORK LAB"/s_series
```

---

## Step 1 — 2D-style structure rendering (single compound)

Render an individual compound as a 3D stick model with ACS-style heteroatom
labels (`NH`, `OH`, `NH₂` rather than bare `N`, `O`).

In PyMOL:

```
load compounds/17sa23.mol, 17sa23
run acs_labels.pml
```

Key points:

- Hydrogen counts are derived from **valence** (N = 3, O = 2, S = 2, minus bond
  orders), not from explicit H atoms, because `.mol` files exported from drawing
  programs usually carry none.
- Aromatic bonds count as order 1.5 so ring nitrogens are typed correctly.
- Labels are sized in Ångströms (`label_size` negative), so they stay sharp at
  any render resolution.
- `float_labels, 1` keeps label backgrounds on top of the sticks instead of
  being clipped by them.

Output: a single PNG of the compound.

---

## Step 2 — Binding-pose figures (all compounds)

```bash
pymol -cq make_docking_figures.py            # focused docking (default)
pymol -cq make_docking_figures.py -- regular # non-focused run
```

For each compound this loads `7L00_AB_protein_FMN.pdb`, overlays the best pose,
selects the FMN cofactor and all residues within 4.0 Å of the ligand, and
renders two views. The camera is set once from the XCJ crystal ligand site and
reused for every compound, so the images align when assembled into a panel.

Output → `figures_focused/`

| File | Contents |
|---|---|
| `<cmpd>_pocket.png` | Ligand (teal) + FMN (gold) in the FabK cartoon |
| `<cmpd>_closeup.png` | Pocket residues labelled, H-bonds dashed with distances |
| `<cmpd>_interactions.txt` | H-bond partners and distances |
| `<cmpd>_session.pse` | Session for manual re-orientation |
| `overlay_all.png` | All five compounds + crystal ligand XCJ in one pocket |
| `overlay_legend.txt` | Colour key for the overlay |

---

## Step 3 — Full interaction analysis (all compounds)

```bash
pymol -cq make_interaction_figures.py            # focused docking (default)
pymol -cq make_interaction_figures.py -- regular # non-focused run
```

Both scripts run headless (`-cq`) and write everything to disk; nothing opens on
screen.

Hydrogens are added to the receptor at runtime for donor/acceptor typing; input
files are not modified. Ligand aromatic rings are found by cycle detection on
the bond graph and kept only if planar (all atoms within 0.25 Å of the mean
plane). Only residues that actually make contact are drawn and labelled, which
keeps the figure readable.

Detection criteria:

| Interaction | Rendered as | Cutoff |
|---|---|---|
| Hydrogen bond | yellow dashes + distance | ≤ 3.5 Å, donor/acceptor |
| Salt bridge | red dashes | ≤ 4.0 Å |
| Halogen bond | magenta dashes | ≤ 3.5 Å, C–X···O/N |
| π–π stacking | green dashes (centroid–centroid) | ≤ 5.5 Å, parallel < 30° or T-shaped > 60° |
| π–cation | orange dashes | ≤ 6.0 Å to Arg/Lys N |
| Hydrophobic | grey dotted | ≤ 4.0 Å C···C, ≥ 2 contacts per residue |

Output → `interactions_focused/`

| File | Contents |
|---|---|
| `<cmpd>_interactions.png` | Contacting residues only, labelled, colour-coded by type |
| `<cmpd>_interactions.tsv` | Every contact: type, residue, chain, atom names, distance |
| `summary_all_compounds.tsv` | All five compounds in one table |
| `legend.txt` | Colour key for the figure caption |
| `<cmpd>_session.pse` | Session for manual re-orientation |

---

## Step 4 — Label sizing

Default label size was too small to read at 2400 × 1800. Adjust in place and
re-render:

```bash
sed -i 's/s("label_size", -[0-9.]*)/s("label_size", -2.0)/' make_interaction_figures.py
pymol -cq make_interaction_figures.py
```

The `[0-9.]*` pattern matches whatever value is currently set, so the command
works no matter how many times the size has already been changed. Values are in
Ångströms; larger magnitude means larger text. `-2.0` was used for the final
figures.

To adjust one figure by hand instead of re-running everything:

```bash
pymol interactions_focused/17sa23_interactions.pse
```
```
set label_size, -2.0
png 17sa23_final.png, width=2400, height=1800, dpi=300, ray=1
```

Other constants worth tuning, at the top of each script:

```python
HB_CUT   = 3.5    # hydrogen-bond cutoff (Å)
PHOB_CUT = 4.0    # hydrophobic C···C cutoff (Å)
WIDTH, HEIGHT, DPI = 2400, 1800, 300
```

---

## Step 5 — Review

```bash
cat interactions_focused/summary_all_compounds.tsv
xdg-open interactions_focused/17sa23_interactions.png
```

Contacts detected in the focused docking run:

| Compound | H-bond | Halogen | π–π | Hydrophobic |
|---|---|---|---|---|
| 5sa23 | 0 | 0 | 1 | 5 |
| 17sa23 | 1 | 2 | 1 | 5 |
| 20sa23 | 2 | 3 | 2 | 4 |
| 22sa23 | 2 | 0 | 1 | 5 |
| 25sa23 | 1 | 5 | 1 | 3 |

20sa23 forms the most extensive polar network; the fluorine-rich compounds
(25sa23, 20sa23) register the most halogen contacts.

---

## Step 6 — Version control

```bash
git init
git branch -M main
git add .
git commit -m "FabK (7L00) structure analysis: scripts, figures and interaction tables"
```

`.gitignore` excludes regenerable and oversized output: `.pse` sessions,
`.pdbqt` files, multi-pose SDF/MOL2 files, MM-GBSA trajectories and the
`figures_*` directories. The `.tsv` interaction tables and `*_best.mol2` poses
are small and are tracked.

Check sizes before pushing (GitHub warns above 50 MB, blocks at 100 MB):

```bash
git count-objects -vH
find . -size +45M -not -path "./.git/*"
```

---

## Issues encountered

| Problem | Cause | Fix |
|---|---|---|
| Labels read `N` instead of `NH` | `.mol` file had no explicit hydrogens, and `h_add` would not add them | Derive H count from valence minus summed bond orders |
| Labels looked like cut-off stickers | White label backgrounds clipped by sticks | `float_labels, 1` and a white outline instead of a filled background |
| Label text invisible against bonds | Label colour matched the bond colour beneath | Black text with a white halo |
| `AttributeError: 'float' object has no attribute 'update'` | Variable named `b` inside `cmd.iterate` collides with PyMOL's reserved B-factor variable | Rename the namespace variables to `p1` / `p2` |
| `sed` label-size edit silently did nothing | Pattern matched the previous literal value, which had already changed | Match any number with `-[0-9.]*` |

---

## Known limitation

The halogen-bond test is distance-only. A genuine halogen bond also requires a
near-linear C–X···O angle (≈140–180°), and aryl C–F is a weak donor at best, so
some reported F···O contacts are close contacts rather than true halogen bonds.
Add an angle filter before quoting these counts.

---

## Notes

- `.pse` sessions are regenerable and are excluded from version control; see the
  top-level `.gitignore`.
- Scripts call `cmd.reinitialize()` on each compound, so anything already loaded
  in an interactive PyMOL session is cleared. Save your work first when running
  these with `run` instead of `pymol -cq`.
