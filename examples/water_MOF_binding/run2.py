#!/Users/Haoyuan/anaconda3/bin/python

import Classes, WaterModelAssign, FFEnergy

newcrd = WaterModelAssign.water('water_MOF841_TraPPEmu3OH.crd','TIP4PEW',[[71,72,73]],geom=False)
FFEnergy.energy(newcrd,['UFF.ff','TraPPEUA.ff','Water.ff'],'','',[],[71,72,73,74],debug=True)  #4-point water models add one site (74 here), make sure to include it in guest atoms!

