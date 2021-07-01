# FFEnergy
Light-weight calculator of binding energy in host-guest complexes using force fields and electrostatics

Quick start:
(1) make a .crd file of your system--very similar to the .xyz format, just add a column of partial charges in the end and remove the first 2 lines
(2) if your system contains water molecules that need to be modeled with standard water models, run WaterModelAssign.py (see examples inside)
(3) run FFEnergy.py (see examples inside)

Available force fields: UFF
Available water models: SPC, SPC/E, TIP3P, TIP4P, TIP4P-Ew
Available special potentials: 12-6-4 LJ

The special potential (.pair file) will override the default 12-6 LJ (but not Coulomb) potential for selected atom pairs.
