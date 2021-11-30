# FFEnergy
Light-weight calculator (with GUI available) of binding energy in host-guest complexes using force fields and electrostatics.

Quick start (GUI):

(1) make a .crd file of your system--very similar to the .xyz format, just add a column of partial charges in the end and remove the first 2 lines--change some atom names if necessary, especially when you're using non-UFF force fields.

(2) open GUI.py, simply go from top to bottom, click corresponding "Load" or "Set" buttons after finishing each part. Visualization is available if you have ase installed.

Quick start (script):

(1) make a .crd file of your system--very similar to the .xyz format, just add a column of partial charges in the end and remove the first 2 lines--change some atom names if necessary, especially when you're using non-UFF force fields.

(2) if your system contains water molecules that need to be modeled with standard water models, run WaterModelAssign.py to generate a new .crd file (see the examples).

(3) run FFEnergy.py (see the examples).

Available force fields: UFF, AMBER-ff99SB, TraPPE-UA.

Available water models: SPC, SPC/E, TIP3P, TIP4P, TIP4P-Ew.

Available special pair potentials: 12-6-4 LJ, Morse.

Only non-bonded interactions, not planning to add bonded terms.

Notes: 

(1) The special potential (.pair file) will override the default 12-6 LJ (but not Coulomb) potential for selected atom pairs.

(2) When reading multiple force fields (.ff files), if same atom names are used in multiple files, the parameters in the last file will override the previous ones.

(3) Cr, Fe and Ce have different 12-6-4 parameters for different oxidation states.
