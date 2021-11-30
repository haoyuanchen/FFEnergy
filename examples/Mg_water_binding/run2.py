#!/Users/Haoyuan/anaconda3/bin/python

import Classes, WaterModelAssign, FFEnergy

newcrd = WaterModelAssign.water('Mg_water.crd','TIP4PEW',[[2,3,4]],geom=False)
FFEnergy.energy(newcrd,['AMBERFF99SB.ff','Water.ff'],'LJ1264.pair','LJ1264',[[1,2]],[2,3,4,5],debug=True)

