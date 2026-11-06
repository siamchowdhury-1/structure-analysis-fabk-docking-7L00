# ============================================================
# Siam Chowdhury (Beck Scholar); [www.siamchowdhury.com]
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Beck College of Sciences and Mathematics, Arkansas State University, AR, USA
# ----------------------------------------------------------------------------
"""
make_docking_figures.py
Publication-style figures for the 7L00 FabK docking results (s_series).

PATHS: all paths are RELATIVE to the current working directory. Developed in
    ~/Documents/Computational_work/WORK LAB/s_series
If you cloned this elsewhere, cd into the folder holding the receptor and the
docking_* folders before running, or edit RECEPTOR / CRYSTAL / COMPOUNDS below.

Run from the s_series folder:
    pymol -cq make_docking_figures.py                 # focused docking (default)
    pymol -cq make_docking_figures.py -- regular      # docking_<cmpd> folders

Makes, in figures_<mode>/ :
    <cmpd>_pocket.png        ligand in the FabK pocket (cartoon + FMN + residues)
    <cmpd>_closeup.png       close-up: H-bonds (dashed, with distances) + labeled residues
    <cmpd>_interactions.txt  list of H-bond partners and distances
    <cmpd>_session.pse       PyMOL session to adjust the view by hand
    overlay_all.png          all compounds + crystal ligand XCJ in one pocket
"""
import os
import sys
from pymol import cmd

# ---------------- settings you may edit ----------------
COMPOUNDS = ["5sa23", "17sa23", "20sa23", "22sa23", "25sa23"]
RECEPTOR = "7L00_AB_protein_FMN.pdb"     # protein + FMN
CRYSTAL = "XCJ_A402_crystal.sdf"          # crystal ligand, sets the common view
POCKET_CUTOFF = 4.0                       # residues within this many A of ligand
HBOND_CUTOFF = 3.6                        # A
WIDTH, HEIGHT, DPI = 2400, 1800, 300

LIG_COLOR = [0.00, 0.60, 0.65]            # teal carbons for the docked ligand
FMN_COLOR = [0.95, 0.75, 0.10]            # gold carbons for FMN
RES_COLOR = [0.55, 0.55, 0.58]            # grey carbons for residues
OVERLAY_COLORS = [[0.00, 0.60, 0.65], [0.85, 0.35, 0.10], [0.45, 0.30, 0.75],
                  [0.20, 0.65, 0.25], [0.85, 0.20, 0.50]]
# --------------------------------------------------------

MODE = "regular" if sys.argv[-1] == "regular" else "focused"
OUT = f"figures_{MODE}"

if not os.path.exists(RECEPTOR):
    sys.exit(f"\nERROR: '{RECEPTOR}' not found in {os.getcwd()}\n"
             "You are running from the wrong folder. cd into the directory that\n"
             "holds the receptor PDB and the docking_* folders, then run again.\n"
             "(Or edit RECEPTOR / CRYSTAL / COMPOUNDS at the top of this script.)\n")

os.makedirs(OUT, exist_ok=True)


def folder(c):
    return f"docking_focused_{c}" if MODE == "focused" else f"docking_{c}"


def style():
    cmd.bg_color("white")
    s = cmd.set
    s("ray_opaque_background", 1)
    s("orthoscopic", 1)
    s("antialias", 2)
    s("ray_trace_mode", 1)
    s("ray_trace_gain", 0.05)
    s("ray_trace_color", "grey20")
    s("ray_shadows", 0)
    s("depth_cue", 0)
    s("specular", 0.3)
    s("shininess", 25)
    s("ambient", 0.35)
    s("direct", 0.6)
    s("valence", 1)
    s("valence_size", 0.07)
    s("stick_radius", 0.17)
    s("cartoon_transparency", 0.15)
    s("cartoon_fancy_helices", 1)
    s("cartoon_side_chain_helper", 1)
    # labels: black bold text, white halo, on top, scaled in Angstrom
    s("label_color", "black")
    s("label_font_id", 7)
    s("label_size", -0.65)
    s("float_labels", 1)
    s("label_outline_color", "white")
    # H-bond dashes
    s("dash_color", "grey20")
    s("dash_radius", 0.06)
    s("dash_gap", 0.3)
    s("dash_length", 0.25)
    s("label_distance_digits", 1)
    cmd.set_color("c_lig", LIG_COLOR)
    cmd.set_color("c_fmn", FMN_COLOR)
    cmd.set_color("c_res", RES_COLOR)


def color_by_element(sel, carbon_color):
    cmd.color(carbon_color, f"({sel}) and e. C")
    cmd.color("tv_blue", f"({sel}) and e. N")
    cmd.color("red", f"({sel}) and e. O")
    cmd.color("yelloworange", f"({sel}) and e. S")
    cmd.color("palecyan", f"({sel}) and e. F")
    cmd.color("green", f"({sel}) and e. Cl")
    cmd.color("white", f"({sel}) and e. H")


def label_residues(sel):
    """ACS-style residue labels, e.g. Tyr147 (on the CA atom)."""
    rows = []
    cmd.iterate(f"({sel}) and name CA", "rows.append((chain, resn, resi))",
                space={"rows": rows})
    chains = {r[0] for r in rows}
    for chain, resn, resi in rows:
        txt = f"{resn.capitalize()}{resi}" + (f" ({chain})" if len(chains) > 1 else "")
        cmd.label(f"({sel}) and chain {chain} and resi {resi} and name CA", repr(txt))


def load_common():
    cmd.reinitialize()
    style()
    cmd.load(RECEPTOR, "rec")
    cmd.remove("rec and solvent")
    if os.path.exists(CRYSTAL):
        cmd.load(CRYSTAL, "xcj")
    cmd.hide("everything")


def write_interactions(c, path):
    pairs = cmd.find_pairs("lig and (donor or acceptor)",
                           "(pocket or fmn) and (donor or acceptor)",
                           mode=1, cutoff=HBOND_CUTOFF)
    lines = [f"{c}  H-bond partners (cutoff {HBOND_CUTOFF} A)",
             "ligand_atom\tpartner\tdistance_A"]
    for (m1, i1), (m2, i2) in pairs:
        a, b = {}, {}
        cmd.iterate(f"{m1} and index {i1}", "p1.update(n=name)", space={"p1": a})
        cmd.iterate(f"{m2} and index {i2}",
                    "p2.update(n=name, r=resn, i=resi, ch=chain)", space={"p2": b})
        d = cmd.get_distance(f"{m1} and index {i1}", f"{m2} and index {i2}")
        lines.append(f"{a['n']}\t{b['r'].capitalize()}{b['i']}{b['ch']}:{b['n']}\t{d:.2f}")
    if len(lines) == 2:
        lines.append("(no H-bonds found)")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return len(lines) - 2


def save_png(path):
    cmd.png(path, width=WIDTH, height=HEIGHT, dpi=DPI, ray=1)
    print("  saved", path)


# ---------------- per-compound figures ----------------
ref_view = None
done = []
for c in COMPOUNDS:
    lig_file = os.path.join(folder(c), f"{c}_best.mol2")
    if not os.path.exists(lig_file):
        print(f"[skip] {c}: {lig_file} not found")
        continue
    print(f"[{c}]")
    load_common()
    cmd.load(lig_file, "lig")
    cmd.remove("lig and h. and neighbor e. C")          # keep only N-H / O-H

    cmd.select("fmn", "byres (rec and resn FMN within 6 of lig)")
    cmd.select("pocket",
               f"byres (rec and polymer.protein within {POCKET_CUTOFF} of lig)")

    # common orientation for all compounds (from the crystal ligand site)
    if ref_view is None:
        cmd.orient("xcj" if "xcj" in cmd.get_names() else "lig")
        cmd.zoom("lig or fmn", 6)
        ref_view = cmd.get_view()
    cmd.set_view(ref_view)

    # ---- 1. pocket overview ----
    cmd.show("cartoon", "rec and polymer.protein")
    cmd.color("grey80", "rec and polymer.protein")
    cmd.show("sticks", "lig or fmn")
    color_by_element("lig", "c_lig")
    color_by_element("fmn", "c_fmn")
    cmd.set("cartoon_transparency", 0.15)
    save_png(os.path.join(OUT, f"{c}_pocket.png"))

    # ---- 2. close-up with H-bonds ----
    cmd.set("cartoon_transparency", 0.6)
    cmd.show("sticks", "pocket and not name N+C+O and not h.")
    color_by_element("pocket", "c_res")
    cmd.distance("hb", "lig", "pocket or fmn", HBOND_CUTOFF, mode=2)
    cmd.hide("labels", "hb")
    cmd.show("labels", "hb")
    cmd.color("grey20", "hb")
    label_residues("pocket")
    cmd.label("fmn and name N5", "'FMN'")
    cmd.zoom("lig", 4)
    save_png(os.path.join(OUT, f"{c}_closeup.png"))

    n = write_interactions(c, os.path.join(OUT, f"{c}_interactions.txt"))
    print(f"  {n} H-bond(s) written")
    cmd.save(os.path.join(OUT, f"{c}_session.pse"))
    done.append((c, lig_file))

# ---------------- overlay of all compounds ----------------
if done:
    print("[overlay]")
    load_common()
    for k, (c, f) in enumerate(done):
        cmd.load(f, c)
        cmd.remove(f"{c} and h.")
        cmd.set_color(f"c_{c}", OVERLAY_COLORS[k % len(OVERLAY_COLORS)])
        color_by_element(c, f"c_{c}")
    lig_all = " or ".join(c for c, _ in done)
    cmd.show("sticks", lig_all)
    cmd.set("stick_radius", 0.13)
    if "xcj" in cmd.get_names():
        cmd.remove("xcj and h.")
        cmd.show("sticks", "xcj")
        cmd.set("stick_radius", 0.09, "xcj")
        cmd.color("grey50", "xcj")
    cmd.show("cartoon", "rec and polymer.protein")
    cmd.color("grey85", "rec and polymer.protein")
    cmd.set("cartoon_transparency", 0.5)
    cmd.show("sticks", f"byres (rec and resn FMN within 6 of ({lig_all}))")
    color_by_element(f"rec and resn FMN", "c_fmn")
    if ref_view:
        cmd.set_view(ref_view)
    # legend text in the console / file
    with open(os.path.join(OUT, "overlay_legend.txt"), "w") as f:
        for k, (c, _) in enumerate(done):
            rgb = OVERLAY_COLORS[k % len(OVERLAY_COLORS)]
            f.write(f"{c}\tRGB {rgb}\n")
        f.write("XCJ (crystal)\tthin grey sticks\nFMN\tgold carbons\n")
    save_png(os.path.join(OUT, "overlay_all.png"))
    cmd.save(os.path.join(OUT, "overlay_session.pse"))

print(f"\nDone. Figures are in {OUT}/")
