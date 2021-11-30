This is an example of calculating the binding energy of a water molecule on a MOF node. The geometry and partial charges are obtained from DFT.

There are two runs here. The "run1.py" will calculate the binding energy using standard UFF force field and the water model you choose. You'll likely get a very positive binding energy of about +200 kJ/mol, which is because the standard UFF LJ parameters cannot properly describe the hydrogen bonding between Zr-OH and water.

The "run2.py" will use the simple fix suggested in J.Phys.Chem.C,2020,124,28469, which is to assign the TraPPE alcohol OH LJ parameters to Zr-OH, to get a more reasonable binding energy. You will get a binding energy of about -70 kJ/mol, which agrees well with the DFT result.

Same calculations can be done using the GUI.
