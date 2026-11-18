"""Protein backbone converter: Cartesian <-> (phi, psi, omega, carbonyl O). Standard library only."""
import math

# Ideal geometry (Engh & Huber, 1991)
L_N_CA, L_CA_C, L_C_N, L_C_O = 1.458, 1.525, 1.329, 1.231
A_N_CA_C, A_CA_C_N, A_C_N_CA, A_CA_C_O = 111.2, 116.2, 121.7, 120.5
PRESETS = {"helix": (-57, -47, 180), "strand": (-139, 135, 180),
           "ppii": (-75, 145, 180), "extended": (180, 180, 180)}

def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def scale(a, s): return [x * s for x in a]
def unit(a): n = math.sqrt(dot(a, a)); return [x / n for x in a]
def dist(a, b): return math.sqrt(dot(sub(a, b), sub(a, b)))

def angle(a, b, c):
    u, v = unit(sub(a, b)), unit(sub(c, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(u, v)))))

def dihedral(p0, p1, p2, p3):
    """IUPAC sign: positive = clockwise looking along p1->p2."""
    b0, b1, b2 = sub(p0, p1), unit(sub(p2, p1)), sub(p3, p2)
    v = sub(b0, scale(b1, dot(b0, b1)))
    w = sub(b2, scale(b1, dot(b2, b1)))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))

def place(a, b, c, length, ang, dih):
    """NeRF: position of atom d bonded to c, given a-b-c."""
    t, p = math.radians(ang), math.radians(dih)
    bc = unit(sub(c, b)); n = unit(cross(sub(b, a), bc)); m = cross(n, bc)
    d = [-length * math.cos(t), length * math.sin(t) * math.cos(p), length * math.sin(t) * math.sin(p)]
    return [c[i] + bc[i]*d[0] + m[i]*d[1] + n[i]*d[2] for i in range(3)]

def build(angles):
    """angles: list of (phi, psi, omega); None allowed only where unused. Returns residue dicts."""
    n = len(angles)
    N = [0.0, 0.0, 0.0]; CA = [L_N_CA, 0.0, 0.0]
    t = math.radians(A_N_CA_C)
    C = [L_N_CA - L_CA_C * math.cos(t), L_CA_C * math.sin(t), 0.0]
    res = []
    for i in range(n):
        phi, psi, om = angles[i]
        if i < n - 1 and (psi is None or om is None): raise ValueError(f"residue {i+1}: psi and omega are required")
        if i > 0 and angles[i][0] is None: raise ValueError(f"residue {i+1}: phi is required")
        O = place(N, CA, C, L_C_O, A_CA_C_O, (psi if psi is not None else 0) + 180)
        res.append({"N": N, "CA": CA, "C": C, "O": O})
        if i < n - 1:
            N = place(res[i]["N"], CA, C, L_C_N, A_CA_C_N, psi)
            CA = place(CA, C, N, L_N_CA, A_C_N_CA, om)
            C = place(C, N, CA, L_CA_C, A_N_CA_C, angles[i + 1][0])
    return res

def to_internal(res):
    out = []
    for i, r in enumerate(res):
        nxt = res[i + 1] if i + 1 < len(res) else None
        linked = nxt is not None and dist(r["C"], nxt["N"]) < 2.0  # else chain break
        out.append({
            "phi": dihedral(res[i-1]["C"], r["N"], r["CA"], r["C"]) if i > 0 and dist(res[i-1]["C"], r["N"]) < 2.0 else None,
            "psi": dihedral(r["N"], r["CA"], r["C"], nxt["N"]) if linked else None,
            "omega": dihedral(r["CA"], r["C"], nxt["N"], nxt["CA"]) if linked else None,
            "C=O": dist(r["C"], r["O"]), "CA-C-O": angle(r["CA"], r["C"], r["O"]),
            "N-CA-C-O": dihedral(r["N"], r["CA"], r["C"], r["O"])})
    return out

def parse_pdb(text):
    res, cur, key = [], None, None
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"): break
        if not ln.startswith("ATOM") or ln[16] not in " A": continue
        name = ln[12:16].strip()
        if name not in ("N", "CA", "C", "O"): continue
        k = (ln[21], ln[22:27])
        if k != key: cur, key = {}, k; res.append(cur)
        cur[name] = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
    bad = [i + 1 for i, r in enumerate(res) if len(r) < 4]
    if bad: raise ValueError(f"residues missing N/CA/C/O: {bad}")
    return res

def to_pdb(res):
    lines, s = [], 1
    for i, r in enumerate(res, 1):
        for name in ("N", "CA", "C", "O"):
            x, y, z = r[name]
            lines.append(f"ATOM  {s:5d} {name:<4s} GLY A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>2s}")
            s += 1
    return "\n".join(lines + ["END"]) + "\n"

def parse_table(text):
    """Lines of 'phi psi omega' (comma/space separated); blank or '-' = missing."""
    rows = []
    for ln in text.strip().splitlines():
        p = [x.strip() for x in ln.replace(",", " ").split()]
        p = (p + ["-"] * 3)[:3]
        rows.append(tuple(None if x in ("-", "") else float(x) for x in p))
    return rows
