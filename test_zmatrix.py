import unittest, backbone as bb, zmatrix as zm
class T(unittest.TestCase):
    def test_pdb_roundtrip(self):
        at = zm.parse_pdb_all(open("inputs/cartesiona_coordinate_Peptide.pdb").read())
        back = zm.zmat_to_cart(zm.parse_zmat(zm.format_zmat(zm.cart_to_zmat(at))))
        self.assertLess(zm.max_distance_error([a[1] for a in at], back), 0.01)
    def test_in_file(self):
        rows = zm.parse_zmat(open("inputs/internal_coordinate.in").read()); x = zm.zmat_to_cart(rows)
        for i, z in enumerate(rows[3:], 3):
            self.assertAlmostEqual(bb.dist(x[i], x[z["na"]-1]), z["r"], 4)
            self.assertAlmostEqual((bb.dihedral(x[i], x[z["na"]-1], x[z["nb"]-1], x[z["nc"]-1]) - z["dih"] + 180) % 360 - 180, 0, 3)
if __name__ == "__main__": unittest.main()
