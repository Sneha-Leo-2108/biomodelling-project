import argparse, sys
import backbone as bb

def fmt(v): return "   --  " if v is None else f"{v:7.2f}"

ap = argparse.ArgumentParser(description="Backbone coordinate converter")
sub = ap.add_subparsers(dest="cmd", required=True)
a = sub.add_parser("to-internal", help="PDB -> phi/psi/omega table"); a.add_argument("pdb")
b = sub.add_parser("to-cartesian", help="angles -> PDB")
b.add_argument("--preset", choices=bb.PRESETS); b.add_argument("--n", type=int, default=10)
b.add_argument("--table", help="file with 'phi psi omega' per line"); b.add_argument("--out", default="output.pdb")
args = ap.parse_args()
try:
    if args.cmd == "to-internal":
        data = bb.to_internal(bb.parse_pdb(open(args.pdb).read()))
        print(" res    phi     psi   omega    C=O  CA-C-O  N-CA-C-O")
        for i, d in enumerate(data, 1):
            print(f"{i:4d} {fmt(d['phi'])} {fmt(d['psi'])} {fmt(d['omega'])} {d['C=O']:6.3f} {d['CA-C-O']:7.2f} {d['N-CA-C-O']:9.2f}")
    else:
        if args.table: angles = bb.parse_table(open(args.table).read())
        elif args.preset: angles = [bb.PRESETS[args.preset]] * args.n
        else: sys.exit("give --preset or --table")
        open(args.out, "w").write(bb.to_pdb(bb.build(angles))); print("wrote", args.out)
except (ValueError, OSError) as e:
    sys.exit(f"error: {e}")
