# Peptide Coordinate Converter
Converts peptide / amino-acid structures between Cartesian coordinates (PDB) and internal coordinates.

- All-atom Z-matrix (bond length, bond angle, dihedral per atom): PDB <-> AMBER-style `.in` file
- Backbone phi, psi, omega (+ carbonyl O): PDB <-> angle table

Files: `backbone.py`, `zmatrix.py` (maths), `main.py` (terminal), `index.html` (website), `inputs/` (given input files), `test_*.py`.

Terminal examples:
    python3 main.py to-zmatrix inputs/cartesiona_coordinate_Peptide.pdb
    python3 main.py from-zmatrix inputs/internal_coordinate.in --resname LYS
    python3 main.py angles inputs/cartesiona_coordinate_Peptide.pdb
    python3 main.py build --preset helix --n 10
Tests: python3 -m unittest discover -s . -p "test_*.py"

Student: YOUR NAME (YOUR ID)
