#!/Users/Haoyuan/anaconda3/bin/python

'''
assign water model parameters to selected water molecules in the crd file
can choose to adjust the geometry of the water or not (geom=True/False)

Haoyuan Chen

Model		R(O-H)		A(H-O-H)	Sigma(O)	Epsilon(O)	q(O)		q(H)		R(O-M)		q(M)
SPC		1.0		109.47		3.166		78.15		-0.82		+0.41		N/A		N/A
SPCE		1.0		109.47		3.166		78.15		-0.8476		+0.4238		N/A		N/A
TIP3P		0.9572		104.52		3.151		76.53		-0.834		+0.417		N/A		N/A
TIP4P		0.9572		104.52		3.154		78.0		0.0		+0.52		0.15		-1.04
TIP4PEW		0.9572		104.52		3.164		81.9		0.0		+0.52422	0.125		-1.04844
'''

from Classes import *
from scipy.spatial.transform import Rotation

def crd_parse(crd_file):
    f = open(crd_file, 'r')
    d = f.readlines()
    f.close()
    atoms = []
    for i in range(len(d)):
        atoms.append(atom(d[i],i+1))
    return atoms

def crd_write(atoms,fname):
    f = open(fname, 'w')
    for a in atoms:
        f.write(a.crdwrite()+'\n')
    f.close()

def geom_adjust(atom_o,atom_h1,atom_h2,water_model):  #keep O fixed, adjust H1's distance, adjust H2's distance and angle
    roh = 0.0
    ahoh = 0.0
    if water_model == 'SPC':
        roh = 1.0
        ahoh = 109.47
    elif water_model == 'SPCE':
        roh = 1.0
        ahoh = 109.47
    elif water_model == 'TIP3P':
        roh = 0.9572
        ahoh = 104.52
    elif water_model == 'TIP4P':
        roh = 0.9572
        ahoh = 104.52
    elif water_model == 'TIP4PEW':
        roh = 0.9572
        ahoh = 104.52
    ahoh = ahoh * np.pi / 180.0
    crd_o = np.array([atom_o.x, atom_o.y, atom_o.z])
    crd_h1 = np.array([atom_h1.x, atom_h1.y, atom_h1.z])
    crd_h2 = np.array([atom_h2.x, atom_h2.y, atom_h2.z])
    vec_oh1 = crd_h1 - crd_o
    vec_oh2 = crd_h2 - crd_o
    r_oh1 = np.linalg.norm(vec_oh1)
    r_oh2 = np.linalg.norm(vec_oh2)
    a_h1oh2 = np.arccos(np.dot(vec_oh1/np.linalg.norm(vec_oh1),vec_oh2/np.linalg.norm(vec_oh2)))
    #print(r_oh1,r_oh2,a_h1oh2)  #for testing
    vec_oh1 = vec_oh1 / r_oh1 * roh
    crd_h1 = crd_o + vec_oh1  #H1 done
    oh1h2_planenorm = np.cross(vec_oh1,vec_oh2)
    rot = Rotation.from_rotvec((ahoh-a_h1oh2) * oh1h2_planenorm / np.linalg.norm(oh1h2_planenorm))
    vec_oh2 = rot.apply(vec_oh2)
    r_oh2 = np.linalg.norm(vec_oh2)
    vec_oh2 = vec_oh2 / r_oh2 * roh
    crd_h2 = crd_o + vec_oh2  #H2 done
    r_oh1 = np.linalg.norm(vec_oh1)
    r_oh2 = np.linalg.norm(vec_oh2)
    a_h1oh2 = np.arccos(np.dot(vec_oh1/np.linalg.norm(vec_oh1),vec_oh2/np.linalg.norm(vec_oh2)))
    #print(r_oh1,r_oh2,a_h1oh2)  #for testing
    atom_h1.x = crd_h1[0]
    atom_h1.y = crd_h1[1]
    atom_h1.z = crd_h1[2]
    atom_h2.x = crd_h2[0]
    atom_h2.y = crd_h2[1]
    atom_h2.z = crd_h2[2]

def add_m(atom_o,atom_h1,atom_h2,water_model):
    #if the water doesn't have perfect geometry for the model, this will put M on the bisector of the H-O-H angle
    crd_o = np.array([atom_o.x, atom_o.y, atom_o.z])
    crd_h1 = np.array([atom_h1.x, atom_h1.y, atom_h1.z])
    crd_h2 = np.array([atom_h2.x, atom_h2.y, atom_h2.z])
    vec_oh1 = crd_h1 - crd_o
    vec_oh2 = crd_h2 - crd_o
    vec_oh1 /= np.linalg.norm(vec_oh1)
    vec_oh2 /= np.linalg.norm(vec_oh2)
    vec_om = (vec_oh1+vec_oh2)/2
    chg = 0.0
    if water_model == 'TIP4P':
        vec_om = vec_om/np.linalg.norm(vec_om)*0.15
        chg = -1.04
    elif water_model == 'TIP4PEW':
        vec_om = vec_om/np.linalg.norm(vec_om)*0.125
        chg = -1.04844
    crd_m = crd_o + vec_om
    crd_line_m = '%s %12.6f %12.6f %12.6f %12.6f'%('Mw'+water_model,crd_m[0],crd_m[1],crd_m[2],chg)
    atom_m = atom(crd_line_m)
    #print(pair(atom_o,atom_m).r)  #checking if O-M distance is right
    return atom_m

def water(crd_file,water_model,waters,geom):  #all atom indices are 1-start, like in FFEnergy.py
    atoms = crd_parse(crd_file)
    for w in waters:  #[[O,H,H],[O,H,H],...]
        atoms[w[0]-1].element = 'Ow'+water_model
        atoms[w[1]-1].element = 'Hw'+water_model
        atoms[w[2]-1].element = 'Hw'+water_model
    if water_model == 'SPC':
        atoms[w[0]-1].charge = -0.82
        atoms[w[1]-1].charge = 0.41
        atoms[w[2]-1].charge = 0.41
    elif water_model == 'SPCE':
        atoms[w[0]-1].charge = -0.8476
        atoms[w[1]-1].charge = 0.4238
        atoms[w[2]-1].charge = 0.4238
    elif water_model == 'TIP3P':
        atoms[w[0]-1].charge = -0.834
        atoms[w[1]-1].charge = 0.417
        atoms[w[2]-1].charge = 0.417
    elif water_model == 'TIP4P':
        atoms[w[0]-1].charge = 0.0
        atoms[w[1]-1].charge = 0.52
        atoms[w[2]-1].charge = 0.52
    elif water_model == 'TIP4PEW':
        atoms[w[0]-1].charge = 0.0
        atoms[w[1]-1].charge = 0.52422
        atoms[w[2]-1].charge = 0.52422
    if geom:
        for w in waters:
            geom_adjust(atoms[w[0]-1],atoms[w[1]-1],atoms[w[2]-1],water_model)
    if (water_model == 'TIP4P') or (water_model == 'TIP4PEW'):
        for w in waters:
            atoms.insert(w[2],add_m(atoms[w[0]-1],atoms[w[1]-1],atoms[w[2]-1],water_model))
        #re-assign all the indices, not necessary though
        for i,a in enumerate(atoms):
            a.idx = i+1
    newcrdfile = ''
    if geom:
        newcrdfile = crd_file.replace('.crd','_'+water_model+'_RIGID.crd')
    else:
        newcrdfile = crd_file.replace('.crd','_'+water_model+'.crd')
    crd_write(atoms,newcrdfile)
    return newcrdfile

#single run example
#water('MOF841_Node_Water_2_Re.crd','SPC',[[71,72,73]],geom=False)
#water('MOF841_Node_Water_2_Re_HtoHw.crd','SPC',[[71,72,73]],geom=False)
#water('MOF841_Node_Water_2_Re_HtoHw.crd','TIP4PEW',[[71,72,73]],geom=False)
#water('MgComplex.crd','TIP4PEW',[[13,14,15]],geom=False)
#water('MgComplex.crd','TIP4PEW',[[13,14,15]],geom=True)

