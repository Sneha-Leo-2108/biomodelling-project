# Backbone Coordinate Converter
Converts protein backbone (N, CA, C, O) between Cartesian coordinates and (phi, psi, omega, carbonyl O).
- `backbone.py` core maths (IUPAC dihedrals, NeRF building, Engh & Huber geometry)
- `main.py` terminal tool: `python3 main.py to-internal file.pdb` / `python3 main.py to-cartesian --preset helix --n 10`
- `index.html` static web app (JS port of the same maths) - deploy by importing this folder in Vercel
- `test_backbone.py` round-trip tests: `python3 test_backbone.py`
