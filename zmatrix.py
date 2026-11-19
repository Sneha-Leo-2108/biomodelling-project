"""All-atom Z-matrix (internal coordinates) <-> Cartesian. Reads AMBER-style .in files and PDB."""
from backbone import dist, angle, dihedral, place, add

def parse_zmat(text):
    rows = []
    for ln in text.splitlines():
        p = ln.split()
        if len(p) < 10 or not p[0].isdigit(): continue
        rows.append(dict(name=p[1], type=p[2], na=int(p[4]), nb=int(p[5]), nc=int(p[6]),
                         r=float(p[7]), ang=float(p[8]), dih=float(p[9])))
    if not rows: raise ValueError("no Z-matrix rows found")
    return rows

def zmat_to_cart(rows):
    xyz = []
    for i, z in enumerate(rows):
        if i == 0: p = [0.0, 0.0, 0.0]
        elif i == 1: p = [z["r"], 0.0, 0.0]
        elif i == 2:
            b, c = xyz[z["nb"] - 1], xyz[z["na"] - 1]
            p = place(add(b, [0, 1, 0]), b, c, z["r"], z["ang"], 0.0)
        else:
            p = place(xyz[z["nc"] - 1], xyz[z["nb"] - 1], xyz[z["na"] - 1], z["r"], z["ang"], z["dih"])
        xyz.append(p)
    return xyz

def element(name):
    return next((ch for ch in name if ch.isalpha()), "X")

def write_pdb(names, xyz, resname="UNK"):
    out = []
    for s, (nm, p) in enumerate(zip(names, xyz), 1):
        shown = f" {nm:<3s}" if len(nm) < 4 else nm
        out.append(f"ATOM  {s:5d} {shown} {resname:<3s} A{1:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00           {element(nm):>2s}")
    return "\n".join(out + ["END"]) + "\n"

def zmat_to_pdb(rows, resname="UNK"):
    keep = [i for i, z in enumerate(rows) if z["name"] != "DUMM"]  # drop dummy atoms
    xyz = zmat_to_cart(rows)
    return write_pdb([rows[i]["name"] for i in keep], [xyz[i] for i in keep], resname)

def parse_pdb_all(text):
    atoms = []
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"): break
        if ln[:6] in ("ATOM  ", "HETATM") and ln[16] in " A":
            atoms.append((ln[12:16].strip(), [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    if not atoms: raise ValueError("no ATOM records found")
    return atoms

def cart_to_zmat(atoms):
    P = [a[1] for a in atoms]; rows = []
    for i, (name, _) in enumerate(atoms):
        na = nb = nc = -1; r = a = d = 0.0
        if i >= 1: na = min(range(i), key=lambda j: dist(P[i], P[j]))
        if i >= 2: nb = min((j for j in range(i) if j != na), key=lambda j: dist(P[na], P[j]))
        if i >= 3:
            cands = sorted((j for j in range(i) if j not in (na, nb)), key=lambda j: dist(P[nb], P[j]))
            nc = next((j for j in cands if 10 < angle(P[na], P[nb], P[j]) < 170), cands[0])
        if i >= 1: r = dist(P[i], P[na])
        if i >= 2: a = angle(P[i], P[na], P[nb])
        if i >= 3: d = dihedral(P[i], P[na], P[nb], P[nc])
        rows.append(dict(name=name, na=na + 1, nb=nb + 1, nc=nc + 1, r=r, ang=a, dih=d))
    return rows

def format_zmat(rows):
    return "\n".join(f"{i:4d}  {z['name']:<5s} {element(z['name']):<4s}  M {max(z['na'],0):4d}{max(z['nb'],0):4d}{max(z['nc'],0):4d} "
                     f"{z['r']:9.3f} {z['ang']:9.3f} {z['dih']:9.3f}  0.000000" for i, z in enumerate(rows, 1)) + "\n"

def max_distance_error(P, Q):
    return max(abs(dist(P[i], P[j]) - dist(Q[i], Q[j])) for i in range(len(P)) for j in range(i))
