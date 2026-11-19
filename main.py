import argparse, sys
import backbone as bb, zmatrix as zm

def fmt(v): return "   --  " if v is None else f"{v:7.2f}"

ap = argparse.ArgumentParser(description="Peptide coordinate converter (Cartesian <-> internal)")
s = ap.add_subparsers(dest="cmd", required=True)
a = s.add_parser("to-zmatrix", help="PDB (all atoms) -> Z-matrix .in file")
a.add_argument("pdb"); a.add_argument("--out", default="output.in")
b = s.add_parser("from-zmatrix", help="Z-matrix .in file -> PDB")
b.add_argument("infile"); b.add_argument("--resname", default="UNK"); b.add_argument("--out", default="output.pdb")
c = s.add_parser("angles", help="PDB -> phi/psi/omega table (backbone)"); c.add_argument("pdb")
d = s.add_parser("build", help="phi/psi/omega -> backbone PDB")
d.add_argument("--preset", choices=bb.PRESETS); d.add_argument("--n", type=int, default=10)
d.add_argument("--table"); d.add_argument("--out", default="backbone.pdb")
args = ap.parse_args()
try:
    if args.cmd == "to-zmatrix":
        rows = zm.cart_to_zmat(zm.parse_pdb_all(open(args.pdb).read()))
        open(args.out, "w").write(zm.format_zmat(rows)); print(f"wrote {args.out} ({len(rows)} atoms)")
    elif args.cmd == "from-zmatrix":
        open(args.out, "w").write(zm.zmat_to_pdb(zm.parse_zmat(open(args.infile).read()), args.resname)); print("wrote", args.out)
    elif args.cmd == "angles":
        print(" res    phi     psi   omega    C=O  CA-C-O  N-CA-C-O")
        for i, x in enumerate(bb.to_internal(bb.parse_pdb(open(args.pdb).read())), 1):
            print(f"{i:4d} {fmt(x['phi'])} {fmt(x['psi'])} {fmt(x['omega'])} {x['C=O']:6.3f} {x['CA-C-O']:7.2f} {x['N-CA-C-O']:9.2f}")
    else:
        if args.table: ang = bb.parse_table(open(args.table).read())
        elif args.preset: ang = [bb.PRESETS[args.preset]] * args.n
        else: sys.exit("give --preset or --table")
        open(args.out, "w").write(bb.to_pdb(bb.build(ang))); print("wrote", args.out)
except (ValueError, OSError) as e:
    sys.exit(f"error: {e}")
