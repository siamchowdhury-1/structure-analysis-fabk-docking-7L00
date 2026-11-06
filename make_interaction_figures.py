# ============================================================
# Siam Chowdhury (Beck Scholar); [www.siamchowdhury.com]
# Computational and Medicinal Chemistry
# Mentor: Mohammad Alam PhD., Professor of Chemistry
# [Dr. Alam's Research Team] www.alamresearch.org
# Beck College of Sciences and Mathematics, Arkansas State University, AR, USA
# ----------------------------------------------------------------------------
"""
make_interaction_figures.py
Full interaction maps for the 7L00 FabK docking results (s_series).

PATHS: all paths are RELATIVE to the current working directory. Developed in
    ~/Documents/Computational_work/WORK LAB/s_series
If you cloned this elsewhere, cd into the folder holding the receptor and the
docking_* folders before running, or edit RECEPTOR / CRYSTAL / COMPOUNDS below.

Run from the s_series folder:
    pymol -cq make_interaction_figures.py                # focused docking (default)
    pymol -cq make_interaction_figures.py -- regular     # docking_<cmpd> folders

Detects, for every compound:
    H-bond            yellow dashes   (donor/acceptor <= 3.5 A)
    Halogen bond      magenta dashes  (C-X...O/N, X = F/Cl/Br/I, <= 3.5 A)
    Salt bridge       red dashes      (charged N ... carboxylate O, <= 4.0 A)
    pi-pi stacking    green dashes    (aromatic ring centroids <= 5.5 A)
    pi-cation         orange dashes   (ring centroid ... Arg/Lys N, <= 6.0 A)
    Hydrophobic       grey dots       (ligand C ... side-chain C, <= 4.0 A)

Writes, in interactions_<mode>/ :
    <cmpd>_interactions.png     figure: only the residues that interact, labeled
    <cmpd>_interactions.tsv     every contact with type, residue, atoms, distance
    <cmpd>_session.pse          PyMOL session to adjust the view by hand
    summary_all_compounds.tsv   one row per compound x residue (for the paper)
    legend.txt                  color key
"""
import os
import sys
import math
from collections import defaultdict
from pymol import cmd
from chempy import cpv

# ---------------- settings ----------------
COMPOUNDS = ["5sa23", "17sa23", "20sa23", "22sa23", "25sa23"]
RECEPTOR = "7L00_AB_protein_FMN.pdb"
WIDTH, HEIGHT, DPI = 2400, 1800, 300

HB_CUT = 3.5          # hydrogen bond
XB_CUT = 3.5          # halogen bond
SB_CUT = 4.0          # salt bridge
PI_CUT = 5.5          # pi-pi centroid distance
PICAT_CUT = 6.0       # pi-cation
PHOB_CUT = 4.0        # hydrophobic C...C
PHOB_MIN_PER_RES = 2  # need this many C...C contacts to call a residue hydrophobic

COLORS = {"hbond": ("yellow", [1.00, 0.85, 0.10]),
          "halogen": ("magenta", [0.85, 0.20, 0.75]),
          "saltbridge": ("red", [0.90, 0.15, 0.15]),
          "pistack": ("green", [0.15, 0.70, 0.25]),
          "pication": ("orange", [1.00, 0.55, 0.05]),
          "hydrophobic": ("grey60", [0.55, 0.55, 0.55])}

HPHOB_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR", "CYS"}
AROM_RES = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
            "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}
CATION_ATOMS = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"]}
ANION_ATOMS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
# -------------------------------------------

MODE = "regular" if sys.argv[-1] == "regular" else "focused"
OUT = f"interactions_{MODE}"

if not os.path.exists(RECEPTOR):
    sys.exit(f"\nERROR: '{RECEPTOR}' not found in {os.getcwd()}\n"
             "You are running from the wrong folder. cd into the directory that\n"
             "holds the receptor PDB and the docking_* folders, then run again.\n"
             "(Or edit RECEPTOR / CRYSTAL / COMPOUNDS at the top of this script.)\n")

os.makedirs(OUT, exist_ok=True)


def folder(c):
    return f"docking_focused_{c}" if MODE == "focused" else f"docking_{c}"


# ---------- geometry helpers ----------
def rings_of(obj, max_size=6):
    """Find 5- and 6-membered rings; return the aromatic (planar) ones."""
    m = cmd.get_model(f"{obj} and not h.")
    idx = [a.index for a in m.atom]
    pos = {a.index: a.coord for a in m.atom}
    adj = defaultdict(set)
    for b in m.bond:
        i, j = m.atom[b.index[0]].index, m.atom[b.index[1]].index
        adj[i].add(j)
        adj[j].add(i)

    found = set()
    for start in idx:
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            for nxt in adj[cur]:
                if nxt == start and len(path) in (5, 6):
                    found.add(tuple(sorted(path)))
                elif nxt not in path and len(path) < max_size:
                    stack.append((nxt, path + [nxt]))

    arom = []
    for ring in found:
        pts = [pos[i] for i in ring]
        cen = [sum(p[k] for p in pts) / len(pts) for k in range(3)]
        n = cpv.cross_product(cpv.sub(pts[1], pts[0]), cpv.sub(pts[2], pts[0]))
        if cpv.length(n) < 1e-6:
            continue
        n = cpv.normalize(n)
        # planarity: every atom close to the mean plane
        if max(abs(cpv.dot_product(cpv.sub(p, cen), n)) for p in pts) < 0.25:
            arom.append({"atoms": list(ring), "center": cen, "normal": n})
    return arom


def res_rings(sel_name):
    """Aromatic rings of PHE/TYR/TRP/HIS side chains in a selection."""
    rows = []
    cmd.iterate(f"({sel_name}) and name CA", "rows.append((chain, resn, resi))",
                space={"rows": rows})
    out = []
    for chain, resn, resi in rows:
        if resn not in AROM_RES:
            continue
        names = "+".join(AROM_RES[resn])
        sel = f"({sel_name}) and chain {chain} and resi {resi} and name {names}"
        co = cmd.get_coords(sel)
        if co is None or len(co) < 5:
            continue
        pts = [list(p) for p in co]
        cen = [sum(p[k] for p in pts) / len(pts) for k in range(3)]
        n = cpv.cross_product(cpv.sub(pts[1], pts[0]), cpv.sub(pts[2], pts[0]))
        if cpv.length(n) < 1e-6:
            continue
        out.append({"key": (chain, resn, resi), "center": cen,
                    "normal": cpv.normalize(n), "sel": sel})
    return out


def atom_info(model, index):
    d = {}
    cmd.iterate(f"{model} and index {index}",
                "d.update(n=name, r=resn, i=resi, ch=chain, e=elem)", space={"d": d})
    return d


def label_text(resn, resi, chain, multi_chain):
    t = f"{resn.capitalize()}{resi}"
    return t + f" ({chain})" if multi_chain else t


# ---------- interaction finders ----------
def find_contacts(cmpd):
    """Return list of dicts: type, reskey, lig_atom, res_atom, dist, coords."""
    out = []

    def add(kind, rk, la, ra, dist, p1, p2):
        out.append({"type": kind, "res": rk, "lig_atom": la, "res_atom": ra,
                    "dist": dist, "p1": p1, "p2": p2})

    # --- H-bonds ---
    for (m1, i1), (m2, i2) in cmd.find_pairs("lig and (donor or acceptor)",
                                             "site and (donor or acceptor)",
                                             mode=1, cutoff=HB_CUT):
        a, b = atom_info(m1, i1), atom_info(m2, i2)
        d = cmd.get_distance(f"{m1} and index {i1}", f"{m2} and index {i2}")
        add("hbond", (b["ch"], b["r"], b["i"]), a["n"], b["n"], d,
            f"{m1} and index {i1}", f"{m2} and index {i2}")

    # --- halogen bonds ---
    hal = cmd.get_model("lig and e. F+Cl+Br+I")
    for at in hal.atom:
        near = cmd.get_model(f"site and e. O+N within {XB_CUT} of "
                             f"(lig and index {at.index})")
        for nb in near.atom:
            d = cpv.distance(at.coord, nb.coord)
            add("halogen", (nb.chain, nb.resn, nb.resi), at.name, nb.name, d,
                f"lig and index {at.index}", f"site and index {nb.index}")

    # --- salt bridges ---
    lig_pos = "lig and e. N and formal_charge > 0"
    lig_neg = "lig and (resn * and name OXT) or (lig and e. O and formal_charge < 0)"
    for resn, names in ANION_ATOMS.items():
        sel = f"site and resn {resn} and name {'+'.join(names)}"
        for (m1, i1), (m2, i2) in cmd.find_pairs(lig_pos, sel, cutoff=SB_CUT):
            a, b = atom_info(m1, i1), atom_info(m2, i2)
            d = cmd.get_distance(f"{m1} and index {i1}", f"{m2} and index {i2}")
            add("saltbridge", (b["ch"], b["r"], b["i"]), a["n"], b["n"], d,
                f"{m1} and index {i1}", f"{m2} and index {i2}")
    for resn, names in CATION_ATOMS.items():
        sel = f"site and resn {resn} and name {'+'.join(names)}"
        for (m1, i1), (m2, i2) in cmd.find_pairs(lig_neg, sel, cutoff=SB_CUT):
            a, b = atom_info(m1, i1), atom_info(m2, i2)
            d = cmd.get_distance(f"{m1} and index {i1}", f"{m2} and index {i2}")
            add("saltbridge", (b["ch"], b["r"], b["i"]), a["n"], b["n"], d,
                f"{m1} and index {i1}", f"{m2} and index {i2}")

    # --- pi-pi stacking and pi-cation ---
    lrings = rings_of("lig")
    rrings = res_rings("site")
    for k, lr in enumerate(lrings):
        pseudo = f"lpi{k}"
        cmd.pseudoatom(pseudo, pos=lr["center"])
        cmd.hide("everything", pseudo)
        for rr in rrings:
            d = cpv.distance(lr["center"], rr["center"])
            if d > PI_CUT:
                continue
            ang = math.degrees(math.acos(min(1.0, abs(
                cpv.dot_product(lr["normal"], rr["normal"])))))
            kind = "pistack"   # parallel (<30 deg) or T-shaped (>60 deg)
            if 30 < ang < 60:
                continue
            cen_name = f"rpi_{rr['key'][1]}{rr['key'][2]}"
            cmd.pseudoatom(cen_name, pos=rr["center"])
            cmd.hide("everything", cen_name)
            add(kind, rr["key"], f"ring{k+1}", "ring", d, pseudo, cen_name)
        # pi-cation
        for resn, names in CATION_ATOMS.items():
            sel = f"site and resn {resn} and name {'+'.join(names)}"
            mm = cmd.get_model(sel)
            for at in mm.atom:
                d = cpv.distance(lr["center"], at.coord)
                if d <= PICAT_CUT:
                    add("pication", (at.chain, at.resn, at.resi),
                        f"ring{k+1}", at.name, d, pseudo,
                        f"site and index {at.index}")

    # --- hydrophobic contacts (best C...C per residue) ---
    best = {}
    for (m1, i1), (m2, i2) in cmd.find_pairs(
            "lig and e. C", "site and e. C and sidechain", cutoff=PHOB_CUT):
        a, b = atom_info(m1, i1), atom_info(m2, i2)
        if b["r"] not in HPHOB_RES:
            continue
        d = cmd.get_distance(f"{m1} and index {i1}", f"{m2} and index {i2}")
        rk = (b["ch"], b["r"], b["i"])
        best.setdefault(rk, {"n": 0, "best": None})
        best[rk]["n"] += 1
        if best[rk]["best"] is None or d < best[rk]["best"][2]:
            best[rk]["best"] = (a["n"], b["n"], d,
                                f"{m1} and index {i1}", f"{m2} and index {i2}")
    for rk, v in best.items():
        if v["n"] >= PHOB_MIN_PER_RES:
            la, ra, d, p1, p2 = v["best"]
            add("hydrophobic", rk, la, ra, d, p1, p2)

    return out


# ---------- rendering ----------
def setup_style():
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
    s("ambient", 0.4)
    s("direct", 0.6)
    s("valence", 1)
    s("valence_size", 0.07)
    s("stick_radius", 0.15)
    s("cartoon_transparency", 0.75)
    s("cartoon_side_chain_helper", 1)
    s("label_color", "black")
    s("label_font_id", 7)
    s("label_size", -0.6)
    s("float_labels", 1)
    s("label_outline_color", "white")
    s("label_distance_digits", 1)
    s("dash_round_ends", 0)
    for name, rgb in COLORS.values():
        pass
    cmd.set_color("c_lig", [0.00, 0.60, 0.65])
    cmd.set_color("c_res", [0.55, 0.55, 0.58])
    cmd.set_color("c_fmn", [0.95, 0.75, 0.10])


def color_elems(sel, carbon):
    cmd.color(carbon, f"({sel}) and e. C")
    cmd.color("tv_blue", f"({sel}) and e. N")
    cmd.color("red", f"({sel}) and e. O")
    cmd.color("yelloworange", f"({sel}) and e. S")
    cmd.color("palecyan", f"({sel}) and e. F")
    cmd.color("green", f"({sel}) and e. Cl")
    cmd.color("brown", f"({sel}) and e. Br")


summary = []
ref_view = None

for c in COMPOUNDS:
    lig_file = os.path.join(folder(c), f"{c}_best.mol2")
    if not os.path.exists(lig_file):
        print(f"[skip] {c}: {lig_file} not found")
        continue
    print(f"[{c}]")
    cmd.reinitialize()
    setup_style()
    cmd.load(RECEPTOR, "rec")
    cmd.remove("rec and solvent")
    cmd.load(lig_file, "lig")
    cmd.remove("lig and h. and neighbor e. C")
    cmd.h_add("rec")                       # needed for donor/acceptor detection
    cmd.select("site", "rec within 6 of lig")
    cmd.hide("everything")

    contacts = find_contacts(c)

    # group by residue, keep best distance per (residue, type)
    bykey = {}
    for ct in contacts:
        k = (ct["res"], ct["type"])
        if k not in bykey or ct["dist"] < bykey[k]["dist"]:
            bykey[k] = ct
    contacts = list(bykey.values())
    reskeys = {ct["res"] for ct in contacts}
    multi = len({r[0] for r in reskeys}) > 1

    # ---- draw ----
    cmd.show("cartoon", "rec and polymer.protein")
    cmd.color("grey90", "rec and polymer.protein")
    cmd.show("sticks", "lig")
    color_elems("lig", "c_lig")

    res_sel_parts = [f"(chain {ch} and resi {ri})" for ch, rn, ri in reskeys]
    if res_sel_parts:
        cmd.select("ires", "rec and (" + " or ".join(res_sel_parts) + ")")
        cmd.show("sticks", "ires and (sidechain or name CA) and not h.")
        color_elems("ires", "c_res")
        if cmd.count_atoms("ires and resn FMN"):
            color_elems("ires and resn FMN", "c_fmn")
        for ch, rn, ri in reskeys:
            sel = f"rec and chain {ch} and resi {ri}"
            anchor = f"({sel}) and name CA"
            if cmd.count_atoms(anchor) == 0:          # FMN / cofactor
                anchor = f"({sel}) and name N5"
                if cmd.count_atoms(anchor) == 0:
                    anchor = f"first ({sel})"
            cmd.label(anchor, repr(label_text(rn, ri, ch, multi)))

    for kind, (pname, rgb) in COLORS.items():
        cmd.set_color(f"ic_{kind}", rgb)
    for n, ct in enumerate(contacts):
        obj = f"{ct['type']}_{n}"
        cmd.distance(obj, ct["p1"], ct["p2"])
        cmd.color(f"ic_{ct['type']}", obj)
        cmd.set("dash_radius", 0.05 if ct["type"] == "hydrophobic" else 0.07, obj)
        cmd.set("dash_gap", 0.45 if ct["type"] == "hydrophobic" else 0.3, obj)
        cmd.set("dash_length", 0.12 if ct["type"] == "hydrophobic" else 0.25, obj)
        if ct["type"] == "hydrophobic":
            cmd.hide("labels", obj)       # distances only on polar contacts

    cmd.orient("lig")
    cmd.zoom("lig", 4.5)
    if ref_view is None:
        ref_view = cmd.get_view()
    cmd.png(os.path.join(OUT, f"{c}_interactions.png"),
            width=WIDTH, height=HEIGHT, dpi=DPI, ray=1)
    cmd.save(os.path.join(OUT, f"{c}_session.pse"))

    # ---- tables ----
    order = ["hbond", "saltbridge", "halogen", "pistack", "pication", "hydrophobic"]
    contacts.sort(key=lambda x: (order.index(x["type"]), x["dist"]))
    with open(os.path.join(OUT, f"{c}_interactions.tsv"), "w") as f:
        f.write("type\tresidue\tchain\tligand_atom\tresidue_atom\tdistance_A\n")
        for ct in contacts:
            ch, rn, ri = ct["res"]
            f.write(f"{ct['type']}\t{rn.capitalize()}{ri}\t{ch}\t"
                    f"{ct['lig_atom']}\t{ct['res_atom']}\t{ct['dist']:.2f}\n")
            summary.append((c, ct["type"], f"{rn.capitalize()}{ri}", ch,
                            ct["lig_atom"], ct["res_atom"], ct["dist"]))
    counts = defaultdict(int)
    for ct in contacts:
        counts[ct["type"]] += 1
    print("  " + ", ".join(f"{k}: {v}" for k, v in counts.items()) or "  none")

# ---- combined summary + legend ----
with open(os.path.join(OUT, "summary_all_compounds.tsv"), "w") as f:
    f.write("compound\ttype\tresidue\tchain\tligand_atom\tresidue_atom\tdistance_A\n")
    for row in summary:
        f.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.2f}"
                          for x in row) + "\n")

with open(os.path.join(OUT, "legend.txt"), "w") as f:
    f.write("Interaction color key\n")
    f.write(f"H-bond\t\tyellow dashes\t<= {HB_CUT} A\n")
    f.write(f"Salt bridge\tred dashes\t<= {SB_CUT} A\n")
    f.write(f"Halogen bond\tmagenta dashes\t<= {XB_CUT} A\n")
    f.write(f"pi-pi stacking\tgreen dashes\t<= {PI_CUT} A (centroid)\n")
    f.write(f"pi-cation\torange dashes\t<= {PICAT_CUT} A\n")
    f.write(f"Hydrophobic\tgrey dots\t<= {PHOB_CUT} A, >= {PHOB_MIN_PER_RES} contacts\n")
    f.write("Ligand carbons teal, residue carbons grey, FMN carbons gold.\n")

print(f"\nDone. Figures and tables are in {OUT}/")
