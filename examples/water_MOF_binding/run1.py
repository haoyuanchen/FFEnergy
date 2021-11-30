#!/Users/Haoyuan/anaconda3/bin/python

import Classes, WaterModelAssign, FFEnergy

newcrd = WaterModelAssign.water('water_MOF841.crd','SPC',[[71,72,73]],geom=False)
FFEnergy.energy(newcrd,['UFF.ff','Water.ff'],'','',[],[71,72,73],debug=False)

