import unittest, backbone as bb
class T(unittest.TestCase):
    def test_roundtrip(self):
        for name, (f, s, o) in bb.PRESETS.items():
            res = bb.build([(f, s, o)] * 8)
            back = bb.to_internal(bb.parse_pdb(bb.to_pdb(res)))
            for i, d in enumerate(back):
                if d["phi"] is not None: self.assertAlmostEqual(abs((d["phi"] - f + 180) % 360 - 180), 0, delta=0.1)
                if d["psi"] is not None: self.assertAlmostEqual(abs((d["psi"] - s + 180) % 360 - 180), 0, delta=0.1)
                if d["omega"] is not None: self.assertAlmostEqual(abs((d["omega"] - o + 180) % 360 - 180), 0, delta=0.1)
    def test_bonds(self):
        r = bb.build([(-57, -47, 180)] * 3)
        self.assertAlmostEqual(bb.dist(r[0]["C"], r[1]["N"]), 1.329, 3)
    def test_mirror_sign(self):
        r = bb.build([(-57, -47, 180)] * 3)
        self.assertLess(bb.to_internal(r)[1]["phi"], 0)
if __name__ == "__main__": unittest.main()
